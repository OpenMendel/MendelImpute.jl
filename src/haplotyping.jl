using NullableArrays

"""
    haplopair!(happair, hapscore, M, N)

Calculate the best pair of haplotypes in `H` for each individual in `X` using
sufficient statistics `M` and `N`.

# Input
* `happair`: storage haplotype pair for each individual.
* `hapmin`: minimum offered by the optimal haplotype pair.
* `M`: `d x d` matrix with entries `M[i, j] = 2dot(H[i, :], H[j, :]) +
    sumabs2(H[i, :]) + sumabs2(H[j, :])`, where `H` is the haplotype matrix
    with haplotypes in rows. Only the upper triangular part of `M` is used.
* `N`: `n x d` matrix `2XH'`, where `X` is the genotype matrix with individuals
    in rows.
"""
function haplopair!(
    happair::Tuple{AbstractVector, AbstractVector},
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
                happair[1][k] = i
                happair[2][k] = j
            end
        end
    end
    return nothing

end

"""
    haplopair(X, H)

Calculate the best pair of haplotypes in `H` for each individual in `X`.

# Input
* `X`: `n x p` genotype matrix. Each row is an individual.
* `H`: `d x p` haplotype matrix. Each row is a haplotype.

# Output
* `happair`: haplotype pair. `X[k, :] ≈ H[happair[1][k], :] + H[happair[2][k], :]`
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
    happair  = zeros(Int, n), zeros(Int, n)
    hapscore = zeros(eltype(N), n)
    haplopair!(X, H, M, N, happair, hapscore)
    return happair, hapscore

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
"""
function haplopair!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    M::AbstractMatrix,
    N::AbstractMatrix,
    happair::Tuple{AbstractVector, AbstractVector},
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
    for I in eachindex(N)
        N[I] *= 2
    end
    haplopair!(happair, hapscore, M, N)
    @inbounds for j in 1:p
        @simd for i in 1:n
            hapscore[i] += abs2(X[i, j])
        end
    end
    return nothing

end

"""
    fillmissing!(X, H, haplopair)

Fill missing genotypes in `X` according to haplotypes.

# Input
* `X`: `n x p` genotype matrix. Each row is an individual.
* `H`: `d x p` haplotype matrix. Each row is a haplotype.
* `happair`: pair of haplotypes. `X[k, :] = H[happair[1][k], :] + H[happair[2][k], :]`.

# Output
* `discrepancy`: sum of squared errors between current values in missing genotypes
    and the imputed genotypes.
"""
function fillmissing!(
    X::NullableMatrix,
    H::AbstractMatrix,
    happair::Tuple{AbstractVector, AbstractVector}
    )

    discrepancy = zero(eltype(X.values))
    for j in 1:size(X, 2), i in 1:size(X, 1)
        if X.isnull[i, j]
            tmp = H[happair[1][i], j] + H[happair[2][i], j]
            discrepancy += abs2(X.values[i, j] - tmp)
            X.values[i, j] = tmp
        end
    end
    return discrepancy

end

"""
    fillgeno!(X, H, happair)

Fill in genotypes according to haplotypes.
"""
function fillgeno!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    happair::Tuple{AbstractVector, AbstractVector}
    )

    @inbounds for j in 1:size(X, 2), i in 1:size(X, 1)
        X[i, j] = H[happair[1][i], j] + H[happair[2][i], j]
    end
    return nothing

end


"""
    haploimpute!(X, H, M, N, happair, hapscore, maxiters=5, tolfun=1e-3)

Haplotying of genotype matrix `X` from the pool of haplotypes `H` and impute
missing genotypes in `X` according to haplotypes.
"""
function haploimpute!(
    X::NullableMatrix,
    H::AbstractMatrix,
    M::AbstractMatrix,
    N::AbstractMatrix,
    happair::Tuple{AbstractVector, AbstractVector},
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
    return nothing

end


function haploimpute!(
    X::NullableMatrix,
    H::AbstractMatrix,
    width::Number = 500,
    maxiters::Int  = 1,
    tolfun::Number = 1e-3,
    verbose::Bool = true
    )

    people, snps = size(X)
    haplotypes = size(H, 1)
    # no need for sliding window
    if nsnps < 3width
        # TODO
    end

    # allocate working arrays
    Xwork = X[:, 1:3width] # NullableMatrix
    Xw1    = view(Xwork.values, :, 1:width)
    Xw23   = view(Xwork.values, :, (width + 1):3width)
    X23    = view(X.values, :, (width + 1):3width)
    Hwork = H[:, 1:3width]
    H1    = view(H, :, 1:width)
    M        = zeros(eltype(H), haplotypes, haplotypes)
    N        = zeros(promote_type(eltype(H), eltype(X.values)), people, haplotypes)
    happair  = zeros(Int, people), zeros(Int, people)
    hapscore = zeros(eltype(N), people

    # phase and impute window 1
    if verbose; println("Imputing SNPs 1:$width"); end
    haploimpute!(Xwork, Hwork, M, N, happair, hapscore, maxiters, tolfun)
    fillgeno!(Xw1, H1, happair)
    copy!(Xw23, X23)

    # phase and impute window 2
    if verbose; println("Imputing SNPs $width+1:2$width"); end
    haploimpute!(Xwork, Hwork, M, N, happair, hapscore, maxiters, tolfun)
    fillgeno!(Xw1, H1, happair)
    copy!(Xw23, X23)


end
