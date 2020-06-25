######## THIS FILE IS THE SAME AS `hpalotype_pair.jl`
######## except it computes a number of top matching haplotype pairs
######## and then screens them each based on the observed entries. 

"""
Records optimal-redundant haplotypes for each window. Currently, only the first 1000
haplotype pairs will be saved to reduce search space for dynamic programming. 

Warning: This function is called in a multithreaded loop. If you modify this function
you must check whether imputation accuracy is affected (when run with >1 threads).
"""
function compute_redundant_haplotypes!(
    redundant_haplotypes::Union{Vector{Vector{Vector{T}}}, Vector{OptimalHaplotypeSet}}, 
    Hunique::CompressedHaplotypes, 
    happairs::Vector{Vector{T}}, 
    window::Int;
    fast_method::Bool = false,
    ) where T <: Tuple{Int32, Int32}
    
    people = length(redundant_haplotypes)
    h1_set = Int32[]
    h2_set = Int32[]

    @inbounds for k in 1:people
        for happair in happairs[k]
            Hi_idx = unique_idx_to_complete_idx(happair[1], window, Hunique)
            Hj_idx = unique_idx_to_complete_idx(happair[2], window, Hunique)

            # loop through all haplotypes and find ones that match either of the optimal haplotypes 
            empty!(h1_set)
            empty!(h2_set)
            for (idx, hap) in enumerate(Hunique.CW_typed[window].hapmap)
                hap == Hi_idx && push!(h1_set, idx)
                hap == Hj_idx && push!(h2_set, idx)
            end

            # save first 1000 haplotype pairs
            for h1 in h1_set, h2 in h2_set
                if length(redundant_haplotypes[k][window]) <= 1000 
                    push!(redundant_haplotypes[k][window], (h1, h2))
                else
                    break
                end
            end
            length(redundant_haplotypes[k][window]) > 1000 && break
        end
    end

    return nothing
end

"""
    haplopair(X, H)

Calculate the best pair of haplotypes in `H` for each individual in `X`. Missing data in `X` 
does not have missing data. Missing data is initialized as 2x alternate allele freq.

# Input
* `X`: `p x n` genotype matrix. Each column is an individual.
* `H`: `p * d` haplotype matrix. Each column is a haplotype.

# Output
* `happair`: optimal haplotype pairs. `X[:, k] ≈ H[:, happair[1][k]] + H[:, happair[2][k]]`.
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
"""
function haplopair_screen(
    X::AbstractMatrix{Union{UInt8, Missing}},
    H::BitMatrix
    )

    Xwork = zeros(Float32, size(X, 1), size(X, 2))
    Hwork = convert(Matrix{Float32}, H)
    initXfloat!(X, Xwork)

    p, n     = size(X)
    d        = size(H, 2)
    M        = zeros(Float32, d, d)
    N        = zeros(Float32, n, d)
    happairs = [Tuple{Int, Int}[] for i in 1:n]
    hapscore = zeros(Float32, n)
    sizehint!.(happairs, 100) # will not save > 100 unique haplotype pairs to conserve memory

    # compute top haplotype pairs for each genotype vector
    t1, t2 = haplopair!(Xwork, Hwork, M, N, happairs, hapscore)
    
    # screen for best haplotype pair based on observed entries
    t3 = @elapsed choose_happair!(X, H, happairs, hapscore)

    return happairs, hapscore, t1, t2, t3
end

"""
    haplopair!(X, H, M, N, happair, hapscore)

Calculate the best pair of haplotypes in `H` for each individual in `X`. Overwite
`M` by `M[i, j] = 2dot(H[:, i], H[:, j]) + sumabs2(H[:, i]) + sumabs2(H[:, j])`,
`N` by `2X'H`, `happair` by optimal haplotype pair, and `hapscore` by
objective value from the optimal haplotype pair.

# Input
* `X`: `p x n` genotype matrix. Each column is an individual.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `M`: overwritten by `M[i, j] = 2dot(H[:, i], H[:, j]) + sumabs2(H[:, i]) +
    sumabs2(H[:, j])`.
* `N`: overwritten by `n x d` matrix `2X'H`.
* `happair`: optimal haplotype pair. `X[:, k] ≈ H[:, happair[k, 1]] + H[:, happair[k, 2]]`.
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
"""
function haplopair!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    M::AbstractMatrix,
    N::AbstractMatrix,
    happairs::Vector{Vector{Tuple{Int, Int}}},
    hapscore::AbstractVector
    )

    p, n, d = size(X, 1), size(X, 2), size(H, 2)

    # assemble M (upper triangular only)
    t1 = @elapsed begin 
        mul!(M, Transpose(H), H)
        for j in 1:d, i in 1:(j - 1) # off-diagonal
            M[i, j] = 2M[i, j] + M[i, i] + M[j, j]
        end
        for j in 1:d # diagonal
            M[j, j] *= 4
        end

        # assemble N
        mul!(N, Transpose(X), H)
        @simd for I in eachindex(N)
            N[I] *= 2
        end
    end

    # computational routine
    t2 = @elapsed haplopair!(happairs, hapscore, M, N)

    # supplement the constant terms in objective
    t2 += @elapsed begin @inbounds for j in 1:n
            @simd for i in 1:p
                hapscore[j] += abs2(X[i, j])
            end
        end
    end

    return t1, t2
end

"""
    haplopair!(happair, hapscore, M, N)

Calculate the best pair of haplotypes pairs in the filtered haplotype panel
for each individual in `X` using sufficient statistics `M` and `N`. 

# Note
The best haplotype pairs are column indices of the filtered haplotype panels. 

# Input
* `happair`: optimal haplotype pair for each individual.
* `hapmin`: minimum offered by the optimal haplotype pair.
* `M`: `d x d` matrix with entries `M[i, j] = 2dot(H[:, i], H[:, j]) +
    sumabs2(H[:, i]) + sumabs2(H[:, j])`, where `H` is the haplotype matrix
    with haplotypes in columns. Only the upper triangular part of `M` is used.
* `N`: `n x d` matrix `2X'H`, where `X` is the genotype matrix with individuals
    in columns.
"""
function haplopair!(
    happairs::Vector{Vector{Tuple{Int, Int}}},
    hapmin::AbstractVector{T},
    M::AbstractMatrix{T},
    N::AbstractMatrix{T},
    tol::T = convert(T, 3)
    ) where T <: Real

    n, d = size(N)
    fill!(hapmin, typemax(T))

    @inbounds for k in 1:d, j in 1:k
        Mjk = M[j, k]
        # loop over individuals
        @simd for i in 1:n
            score = Mjk - N[i, j] - N[i, k]

            # keep best happair (original code)
            # if score < hapmin[i]
            #     hapmin[i], happair1[i], happair2[i] = score, j, k
            # end

            # different data structure but equivalent to original code
            # if score < hapmin[i]
            #     empty!(happairs[i])
            #     push!(happairs[i], (j, k))
            #     hapmin[i] = score
            # end

            # keep all happairs that are equally good
            # if score < hapmin[i]
            #     empty!(happairs[i])
            #     push!(happairs[i], (j, k))
            #     hapmin[i] = score
            # elseif score == hapmin[i]
            #     push!(happairs[i], (j, k))
            # end

            # keep happairs that within some range of best pair (but finding all of them requires a 2nd pass)
            if score < hapmin[i]
                empty!(happairs[i])
                push!(happairs[i], (j, k))
                hapmin[i] = score
            elseif score <= hapmin[i] + tol && length(happairs[i]) < 100
                push!(happairs[i], (j, k))
            end

            # keep top 10 haplotype pairs
            # if score < hapmin[i]
            #     length(happairs[i]) == 10 && popfirst!(happairs[i])
            #     push!(happairs[i], (j, k))
            #     hapmin[i] = score
            # elseif score <= hapmin[i] + interval
            #     length(happairs[i]) == 10 && popfirst!(happairs[i])
            #     push!(happairs[i], (j, k))
            # end

            # keep all previous best pairs
            # if score < hapmin[i]
            #     push!(happairs[i], (j, k))
            #     hapmin[i] = score
            # end

            # keep all previous best pairs and equally good pairs
            # if score <= hapmin[i]
            #     push!(happairs[i], (j, k))
            #     hapmin[i] = score
            # end
        end
    end

    return nothing
end

"""
    haploimpute!(X, H, M, N, happair, hapscore, maxiters=1, tolfun=1e-3)

In a window, performs haplotying of genotype matrix `X` from the pool of 
haplotypes `H`.

# Input
* `X`: `p x n` matrix with missing values. Each column is genotypes of an individual.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `M`: overwritten by `M[i, j] = 2dot(H[:, i], H[:, j]) + sumabs2(H[:, i]) +
    sumabs2(H[:, j])`.
* `N`: overwritten by `n x d` matrix `2X'H`.
* `happair`: vector of optimal haplotype pair. `X[:, k] ≈ H[:, happair[k, 1]] + H[:, happair[k, 2]]`.
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
* `Xfloat`: copy of `X` where missing values are filled with mean. This engages in linear algebra for computing `N`
* `maxiters`: number of MM iterations. Default is 1.
* `tolfun`: convergence tolerance of MM iterations. Default is 1e-3.
"""
function haploimpute!(
    X::AbstractMatrix,
    H::AbstractMatrix{T},
    M::AbstractMatrix{T},
    N::AbstractMatrix{T},
    happairs::Vector{Vector{Tuple{Int, Int}}},
    hapscore::AbstractVector;
    Xfloat::AbstractMatrix = zeros(T, size(X)),
    maxiters::Int  = 1,
    tolfun::Number = 1e-3,
    ) where T

    obj = typemax(eltype(hapscore))
    size(X) == size(Xfloat) || error("Dimension mismatch: X and Xfloat have sizes $(size(X)) and $(size(Xfloat))")
    initXfloat!(X, Xfloat) #Xfloat is the matrix that engages in BLAS routines

    # mm iteration
    for iter in 1:maxiters
        # compute top haplotype pairs for each genotype vector
        haplopair!(Xfloat, H, M, N, happairs, hapscore)
        # screen for best haplotype pair based on observed entries
        choose_happair!(X, H, happairs, hapscore)
        # impute missing entries according to current haplotypes
        discrepancy = fillmissing!(X, Xfloat, H, happairs)
        # convergence criterion
        objold = obj
        obj = sum(hapscore) - discrepancy
        # println("iter = $iter, discrepancy = $discrepancy, obj = $obj")
        if abs(obj - objold) < tolfun * (objold + 1)
            break
        end
    end
    return nothing
end

"""
    choose_happair!(X, H, happairs, hapscore)

Calculates error ||x - hi - hj||^2 only on the observed entries and save 
observed error in `hapscore`. `happairs` will keep only the best haplotype 
pairs based on the error of observed entries. All happairs
that attain the best observed error will be kept.
"""
function choose_happair!(
    X::AbstractMatrix{Union{Missing, T}},
    H::AbstractMatrix,
    happairs::Vector{Vector{Tuple{Int, Int}}},
    hapscore::AbstractVector;
    ) where T <: Real

    p = size(X, 1)
    n = size(X, 2)
    d = size(H, 2)
    p == size(H, 1) || error("Dimension mismatch: size(X, 1) = $p but size(H, 1) = $(size(H, 1))")

    # loop over each person's genotype
    best_happair = Tuple{Int, Int}[]
    for j in 1:n
        best_error = typemax(T)
        empty!(best_happair)
        for happair in happairs[j]
            # compute errors for each pair based on observed entries
            h1, h2 = happair[1], happair[2]
            err = zero(T)
            @inbounds @simd for i in 1:p
                if X[i, j] !== missing 
                    err += abs2(X[i, j] - H[i, h1] - H[i, h2])
                end
            end
            if err :: T < best_error
                best_error = err
                empty!(best_happair)
                push!(best_happair, happair)
            elseif err :: T == best_error
                push!(best_happair, happair)
            end
        end

        # keep only best haplotype pair in happairs
        if length(happairs[j]) > 1
            empty!(happairs[j])
            for pair in best_happair
                push!(happairs[j], pair)
            end
        end
        hapscore[j] = convert(eltype(hapscore), best_error)
    end

    return nothing
end
