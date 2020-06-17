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
    ) where T <: Tuple{Int, Int}
    
    people = length(redundant_haplotypes)
    h1_set = Int[]
    h2_set = Int[]

    @inbounds for k in 1:people, happair in happairs[k]
        Hi_idx = unique_idx_to_complete_idx(happair[1], window, Hunique)
        Hj_idx = unique_idx_to_complete_idx(happair[2], window, Hunique)

        # loop through all haplotypes and find ones that match either of the optimal haplotypes 
        empty!(h1_set)
        empty!(h2_set)
        for (idx, hap) in enumerate(Hunique[window].hapmap)
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
function haplopair(
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
    haplopair!(Xwork, Hwork, M, N, happairs, hapscore)

    return happairs, hapscore
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

    # computational routine
    haplopair!(happairs, hapscore, M, N)

    # supplement the constant terms in objective
    @inbounds for j in 1:n
        @simd for i in 1:p
            hapscore[j] += abs2(X[i, j])
        end
    end

    return nothing
end

"""
    haplopair!(happair, hapscore, M, N)

Calculate the best pair of haplotypes pairs in the filtered haplotype panel
for each individual in `X` using sufficient statistics `M` and `N`. 

# Note
The best haplotype pairs are column indices of the filtered haplotype 
panel, and must be converted back to the indices of the actual haplotype panel
using `UniqueHaplotypeMaps`. 

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
    hapmin::Vector{T},
    M::AbstractMatrix{T},
    N::AbstractMatrix{T},
    ) where T <: Real

    n, d = size(N)
    tol = convert(T, 3)
    fill!(hapmin, typemax(T))
    empty!.(happairs)

    @inbounds for k in 1:d, j in 1:k
        # loop over individuals
        @simd for i in 1:n
            score = M[j, k] - N[i, j] - N[i, k]

            # keep best happair (original code)
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
    fillmissing!(Xm, Xwork, H, haplopairs)

Fill in missing genotypes in `X` according to haplotypes. Non-missing genotypes
remain same.

# Input
* `Xm`: `p x n` genotype matrix with missing values. Each column is an individual.
* `Xwork`: `p x n` genotype matrix where missing values are filled with sum of 2 haplotypes.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `happair`: pair of haplotypes. `X[:, k] = H[:, happair[1][k]] + H[:, happair[2][k]]`.
"""
function fillmissing!(
    Xm::AbstractMatrix{Union{U, Missing}},
    Xwork::AbstractMatrix{T},
    H::AbstractMatrix{T},
    happairs::Vector{Vector{Tuple{Int, Int}}},
    ) where {T, U}

    p, n = size(Xm)
    best_discrepancy = typemax(eltype(Xwork))
    
    for j in 1:n, happair in happairs[j]
        discrepancy = zero(T)
        for i in 1:p
            if ismissing(Xm[i, j])
                tmp = H[i, happair[1]] + H[i, happair[2]]
                discrepancy += abs2(Xwork[i, j] - tmp)
                Xwork[i, j] = tmp
            end
        end
        if discrepancy < best_discrepancy
            best_discrepancy = discrepancy
        end
    end
    return best_discrepancy
end

# """
#     fillgeno!(X, H, happair)

# Fill in genotypes according to haplotypes. Both missing and non-missing
# genotypes may be changed.

# # Input
# * `X`: `p x n` genotype matrix. Each column is an individual.
# * `H`: `p x d` haplotype matrix. Each column is a haplotype.
# * `happair`: pair of haplotypes. `X[:, k] = H[:, happair[1][k]] + H[:, happair[2][k]]`.
# """
# function fillgeno!(
#     X::AbstractMatrix,
#     H::AbstractMatrix,
#     happair::Tuple{AbstractVector, AbstractVector}
#     )

#     @inbounds for j in 1:size(X, 2), i in 1:size(X, 1)
#         X[i, j] = H[i, happair[1][j]] + H[i, happair[2][j]]
#     end
#     return nothing

# end

"""
    initXfloat!(X, Xfloat)

Initializes the matrix `Xfloat` where missing values of matrix `X` by `2 x` allele frequency
and nonmissing entries of `X` are converted to type `Float32` for subsequent BLAS routines. 

# Input
* `X` is a `p x n` genotype matrix. Each column is an individual.
* `Xfloat` is the `p x n` matrix of X where missing values are filled by 2x allele frequency. 
"""
function initXfloat!(
    X::AbstractMatrix,
    Xfloat::AbstractMatrix
    )
    
    T = Float32
    p, n = size(X)

    @inbounds for i in 1:p
        # allele frequency
        cnnz = zero(T)
        csum = zero(T)
        for j in 1:n
            if !ismissing(X[i, j])
                cnnz += one(T)
                csum += convert(T, X[i, j])
            end
        end
        # set missing values to 2freq, unless cnnz is 0
        imp = (cnnz == 0 ? zero(T) : csum / cnnz)
        for j in 1:n
            if ismissing(X[i, j])
                Xfloat[i, j] = imp
            else
                Xfloat[i, j] = convert(T, X[i, j])
            end
        end
    end

    any(isnan, Xfloat) && error("Xfloat contains NaN during initialization! Shouldn't happen!")
    any(isinf, Xfloat) && error("Xfloat contains Inf during initialization! Shouldn't happen!")
    any(ismissing, Xfloat) && error("Xfloat contains Missing during initialization! Shouldn't happen!")

    # impute using mode
    # for i in 1:p
    #     # set missing values to mode
    #     imp = mode(@view(X[i, :]))
    #     for j in 1:n
    #         if ismissing(X[i, j]) 
    #             Xfloat[i, j] = imp
    #         else
    #             Xfloat[i, j] = X[i, j]
    #         end
    #     end
    # end

    # initialize using 0
    # for i in 1:p, j in 1:n
    #     Xfloat[i, j] = ifelse(ismissing(X[i, j]), zero(T), X[i, j])
    # end

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
    U = promote_type(eltype(X), eltype(H))

    # loop over each person's genotype
    best_happair = Tuple{Int, Int}[]
    for j in 1:n
        best_error = typemax(eltype(H))
        empty!(best_happair)
        for happair in happairs[j]
            # compute errors for each pair based on observed entries
            h1, h2 = happair[1], happair[2]
            err = zero(U)
            @inbounds @simd for i in 1:p
                if X[i, j] !== missing 
                    err += abs2(X[i, j] - H[i, h1] - H[i, h2])
                end
            end
            if err :: U < best_error
                best_error = err
                empty!(best_happair)
                push!(best_happair, happair)
            elseif err :: U == best_error
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

"""
    chunk_size(people, haplotypes)

Figures out how many SNPs can be loaded into memory (capped at 2/3 total RAM) 
at once, given the data size. Assumes genotype data are Float32 (4 byte per entry) 
and haplotype panels are BitArrays (1 bit per entry).
"""
function chunk_size(people::Int, haplotypes::Int)
    system_memory_gb = Sys.total_memory() / 2^30
    system_memory_bits = 8000000000 * system_memory_gb
    usable_bits = round(Int, system_memory_bits * 2 / 3) # use 2/3 of memory for genotype and haplotype matrix per chunk
    max_chunk_size = round(Int, usable_bits / (haplotypes + 32people))
    return max_chunk_size

    # return 1000 # for testing in compare1
end
