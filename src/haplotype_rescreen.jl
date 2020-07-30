######## THIS FILE IS THE SAME AS `haplotype_pair.jl`
######## except it computes a number of top matching haplotype pairs
######## and then screens them each based on the observed entries.

"""
    haplopair_rescreen(X, H, ...)

Calculate the best pair of haplotypes in `H` for each individual in `X`.
Missing data in `X` does not have missing data. Missing data is initialized as
2x alternate allele freq.

# Input
* `X`: `p x n` genotype matrix. Each column is an individual.
* `H`: `p * d` haplotype matrix. Each column is a haplotype.

# Output
* `happair`: optimal haplotype pairs. `X[:, k] ≈ H[:, happair[1][k]] + 
    H[:, happair[2][k]]`.
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
"""
function haplopair_rescreen!(
    X::AbstractMatrix,
    H::BitMatrix;
    # preallocated vectors
    happair1::AbstractVector = ones(Int, size(X, 2)),      # length n
    happair2::AbstractVector = ones(Int, size(X, 2)),      # length n
    hapscore::AbstractVector = Vector{Float32}(undef, size(X, 2)), # length n
    # preallocated matrices
    M     :: AbstractMatrix{Float32} = Matrix{Float32}(undef, size(H, 2), size(H, 2)), # cannot be preallocated until Julia 2.0
    Xwork :: AbstractMatrix{Float32} = Matrix{Float32}(undef, size(X, 1), size(X, 2)), # p × n
    Hwork :: AbstractMatrix{Float32} = convert(Matrix{Float32}, H),                    # p × d (not preallocated)
    N     :: AbstractMatrix{Float32} = Matrix{Float32}(undef, size(X, 2), size(H, 2)), # n × d (not preallocated)
    )

    p, n = size(X)
    d    = size(H, 2)

    # reallocate matrices for last window
    if size(Xwork, 1) != p
        Xwork = zeros(Float32, p, n)
    end

    # initialize missings in Xwork
    initXfloat!(Xwork, X)

    # working array
    happairs = [Tuple{Int32, Int32}[] for i in 1:n]
    sizehint!.(happairs, 100) # save only first 100 pairs to conserve memory

    # compute top haplotype pairs for each genotype vector
    t2, t3 = haplopair!(Xwork, Hwork, M, N, happairs, hapscore)

    # screen for best haplotype pair based on observed entries
    t4 = @elapsed choose_happair!(X, H, happairs, hapscore)

    # record best haplotype based on observed entries
    for i in 1:n
        best = happairs[i][1]
        happair1[i] = best[1]
        happair2[i] = best[2]
    end

    t1 = 0.0 # no time spent on haplotype thinning
    return t1, t2, t3, t4
end

"""
    haplopair!(X, H, M, N, happair, hapscore)

Calculate the best pair of haplotypes in `H` for each individual in `X`.
Overwite `M` by `M[i, j] = 2dot(H[:, i], H[:, j]) + sumabs2(H[:, i]) + 
sumabs2(H[:, j])`, `N` by `2X'H`, `happair` by optimal haplotype pair, and 
`hapscore` by objective value from the optimal haplotype pair.

# Input
* `X`: `p x n` genotype matrix. Each column is an individual.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `M`: overwritten by `M[i, j] = 2dot(H[:, i], H[:, j]) + sumabs2(H[:, i]) +
    sumabs2(H[:, j])`.
* `N`: overwritten by `n x d` matrix `2X'H`.
* `happair`: optimal haplotype pair. `X[:, k] ≈ H[:, happair[k, 1]] + 
    H[:, happair[k, 2]]`.
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
"""
function haplopair!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    M::AbstractMatrix,
    N::AbstractMatrix,
    happairs::Vector{Vector{Tuple{Int32, Int32}}},
    hapscore::AbstractVector
    )

    p, n, d = size(X, 1), size(X, 2), size(H, 2)

    # assemble M (upper triangular only)
    t2 = @elapsed begin
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
    t3 = @elapsed haplopair!(happairs, hapscore, M, N)

    # supplement the constant terms in objective
    t3 += @elapsed begin @inbounds for j in 1:n
            @simd for i in 1:p
                hapscore[j] += abs2(X[i, j])
            end
        end
    end

    return t2, t3
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
    happairs::Vector{Vector{Tuple{Int32, Int32}}},
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

            # keep happairs that within some range of best pair
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
Computes ||X[:, col] - H[:, h1] - H[:, h2] ||^2 on just the observed entries.
"""
function observed_error(X, col, H, h1, h2)
    p = size(X, 1)
    @assert p == size(H, 1)
    T = promote_type(eltype(X), eltype(H))
    err = zero(T)
    @inbounds for i in 1:p
        if X[i, col] !== missing
            err += abs2(X[i, col] - H[i, h1] - H[i, h2])
        end
    end
    return err :: T
end

"""
    choose_happair!(X, H, happairs, hapscore)

Calculates error ||x - hi - hj||^2 only on the observed entries and save
observed error in `hapscore`. `happairs` will keep only the best haplotype
pairs based on the error of observed entries. All happairs
that attain the best observed error will be kept.
"""
function choose_happair!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    happairs::Vector{Vector{Tuple{Int32, Int32}}},
    hapscore::AbstractVector;
    )

    p = size(X, 1)
    n = size(X, 2)
    d = size(H, 2)
    p == size(H, 1) || error("Dimension mismatch: size(X, 1) = $p but" *
        " size(H, 1) = $(size(H, 1))")
    T = UInt8 # X is Matrix{Union{Missing, UInt8}}

    # loop over each person's genotype
    best_h1 = best_h2 = 0
    for j in 1:n
        # compute errors for each pair based on observed entries
        best_error = typemax(T)
        for happair in happairs[j]
            h1, h2 = happair[1], happair[2]
            err = observed_error(X, j, H, h1, h2)
            if err < best_error
                best_error, best_h1, best_h2 = err, h1, h2
            end
        end

        # keep only best haplotype pair in happairs
        empty!(happairs[j])
        push!(happairs[j], (best_h1, best_h2))
        hapscore[j] = convert(eltype(hapscore), best_error)
    end

    return nothing
end
