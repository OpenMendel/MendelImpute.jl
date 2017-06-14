"""
    haplopair!(happair, hapscore, M, N)

Calculate the best pair of haplotypes in `H` for each individual in `X` using
sufficient statistics `M` and `N`.

# Input
* `happair`: `n * 2` storage for haplotype index for each individual.
* `hapmin`: minimum offered by the optimal haplotype pair.
* `M`: `d x d` matrix with entries `M[i, j] = dot(H[i, :], H[j, :]) +
    sumabs2(H[i, :]) + sumabs2(H[j, :])`, where `H` is the haplotype matrix
    with haplotypes in rows. Only the upper triangular part of `M` is used.
* `N`: `n x d` matrix `2XH'`, where `X` is the genotype matrix with individuals
    in rows.

# Output
* `happair`
* `hapmin`
"""
function haplopair!(
    happair::AbstractMatrix,
    hapmin::Vector,
    M::AbstractMatrix,
    N::AbstractMatrix
    )

    n, d = size(N)
    fill!(hapmin, typemax(eltype(hapmin)))
    @inbounds for j in 1:d, i in 1:j
        mij = M[i, j]
        # loop over individuals
        @simd for k in 1:n
            score = mij - N[k, i] - N[k, j]
            if score < hapmin[k]
                hapmin[k]     = score
                happair[k, 1] = i
                happair[k, 2] = j
            end
        end
    end
    return happair, hapmin

end

"""
    haplopair(X, H)

Calculate the best pair of haplotypes in `H` for each individual in `X`.

# Input
* `X`: `n x p` genotype matrix. Each row is an individual.
* `H`: `d x p` haplotype matrix. Each row is a haplotype.

# Output
* `happair`: haplotype pair. `X[k, :] ≈ H[happair[k, 1], :] + H[happair[k, 2], :]`
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
"""
function haplopair(
    X::AbstractMatrix,
    H::AbstractMatrix
    )

    n, p     = size(X)
    d        = size(H, 1)
    M        = zeros(eltype(H), d, d)
    N        = zeros(promote_type(eltype(H), eltype(X)), n, d)
    happair  = zeros(Int, n, 2)
    hapscore = zeros(eltype(N), n)
    haplopair!(X, H, M, N, happair, hapscore)

end

"""
    haplopair!(X, H, M, N, happair, hapscore)

Calculate the best pair of haplotypes in `H` for each individual in `X`. Overwite
`M` by `2HH'`, `N` by `2XH'`, `happair` by optimal haplotype pair, and `hapscore`
by objective value from optimal haplotype pair.

# Input
* `X`: `n x p` genotype matrix. Each row is an individual.
* `H`: `d x p` haplotype matrix. Each row is a haplotype.
* `M`: overwritten by `d x d` matrix `2HH'` with diagonal `sumabs2(H[i, :])`.
* `N`: overwritten by `n x d` matrix `2XH'`.
* `happair`: haplotype pair. `X[k, :] ≈ H[happair[k, 1], :] + H[happair[k, 2], :]`
* `hapscore`: haplotyping score. 0 means best. Larger means worse.

# Output
* `happair`
* `hapscore`
"""
function haplopair!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    M::AbstractMatrix,
    N::AbstractMatrix,
    happair::AbstractMatrix,
    hapscore::AbstractVector
    )

    n, p = size(X)
    d    = size(H, 1)
    A_mul_Bt!(M, H, H)
    for j in 1:d
        for i in 1:(j - 1)
            M[i, j] = 2M[i, j] + M[i, i] + M[j, j]
        end
    end
    for j in 1:d
        M[j, j] *= 4
    end
    A_mul_Bt!(N, X, H)
    N .*= 2
    haplopair!(happair, hapscore, M, N)
    @inbounds for j in 1:p
        @simd for i in 1:n
            hapscore[i] += X[i, j] * X[i, j]
        end
    end
    return happair, hapscore

end

"""
    fillmissing!(X, H, haplopair)

Fill missing genotypes in `X` according to haplotypes.

# Input
* `X`: `n x p` genotype matrix. Each row is an individual.
* `H`: `d x p` haplotype matrix. Each row is a haplotype.
* `happair`: pair of haplotypes. `X[k, :] = H[happair[k, 1], :] + H[happair[k, 2], :]`.
* `missingidx`: `X[missingidx[1][i], missingidx[2][i]]` are missing.

# Output
* `X`: completed genotype matrix.
"""
function fillmissing!(
    X::NullableMatrix,
    H::AbstractMatrix,
    happair::AbstractMatrix,
    )

    n, p = size(X)
    discrepancy = zero(eltype(X))
    for j in 1:size(X, 2)
        for i in 1:size(X, 1)
            if X.isnull[i, j]
                tmp = H[happair[i, 1], j] + H[happair[i, 2], j]
                discrepancy += (X.values[i, j] - tmp)^2
                X.values[i, j] = tmp
            end
        end
    end
    return discrepancy

end

"""
    haploimpute!(X, H, M, N, happair, hapscore, maxiters=5, tolfun=1e-3)
"""
function haploimpute!(
    X::NullableMatrix,
    H::AbstractMatrix,
    M::AbstractMatrix,
    N::AbstractMatrix,
    happair::AbstractMatrix,
    hapscore::AbstractVector,
    maxiters::Int  = 1,
    tolfun::Number = 1e-3
    )

    obj = typemax(eltype(hapscore))
    for iter in 1:maxiters
        # haplotyping
        haplopair!(X.values, H, M, N, happair, hapscore)
        # impute missing entries according to current haplotypes
        discrepancy = fillmissing!(X, H, happair)
        #println("discrepancy = $discrepancy")
        # convergence criterion
        objold = obj
        obj = sum(hapscore) - discrepancy
        #println("iter = $iter, obj = $obj")
        if abs(obj - objold) < tolfun * (objold + 1)
            break
        end
    end

end
