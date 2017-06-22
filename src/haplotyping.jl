using NullableArrays
import Base.Order: Lexicographic

"""
    haplopair!(happair, hapscore, M, N)

Calculate the best pair of haplotypes in `H` for each individual in `X` using
sufficient statistics `M` and `N`.

# Input
* `happair`: optimal haplotype pair for each individual.
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
    # TODO: parallel computing
    @inbounds for j in 1:d, i in 1:j
        mij = M[i, j]
        # loop over individuals
        for k in 1:n
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
`M` by `M[i, j] = 2dot(H[i, :], H[j, :]) + sumabs2(H[i, :]) + sumabs2(H[j, :])`,
`N` by `2XH'`, `happair` by optimal haplotype pair, and `hapscore` by
objective value from the optimal haplotype pair.

# Input
* `X`: `n x p` genotype matrix. Each row is an individual.
* `H`: `d x p` haplotype matrix. Each row is a haplotype.
* `M`: overwritten by `M[i, j] = 2dot(H[i, :], H[j, :]) + sumabs2(H[i, :]) +
    sumabs2(H[j, :])`.
* `N`: overwritten by `n x d` matrix `2XH'`.
* `happair`: optimal haplotype pair. `X[k, :] ≈ H[happair[k, 1], :] + H[happair[k, 2], :]`
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
    for j in 1:d, i in 1:(j - 1)
        M[i, j] = 2M[i, j] + M[i, i] + M[j, j]
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

Fill in missing genotypes in `X` according to haplotypes. Non-missing gentypes
remain same.

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

    discrepancy = zero(promote_type(eltype(X.values), eltype(H)))
    @inbounds for j in 1:size(X, 2), i in 1:size(X, 1)
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

Fill in genotypes according to haplotypes. Both missing and non-missing
genotypes may be changed.

# Input
* `X`: `n x p` genotype matrix. Each row is an individual.
* `H`: `d x p` haplotype matrix. Each row is a haplotype.
* `happair`: pair of haplotypes. `X[k, :] = H[happair[1][k], :] + H[happair[2][k], :]`.
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
    initmissing(X)

Initialize the missing values in a nullable matrix `X` by `2 x` allele frequency.
"""
function initmissing!(X::NullableMatrix)
    T = eltype(X.values)
    for j in 1:size(X, 2)
        # allele frequency
        cnnz = 0
        csum = zero(T)
        for i in 1:size(X, 1)
            if ~X.isnull[i, j]
                cnnz += 1
                csum += X.values[i, j]
            end
        end
        # set missing values to 2freq
        imp = csum / cnnz
        for i in 1:size(X, 1)
            if X.isnull[i, j]
                X.values[i, j] = imp
            end
        end
    end
    return nothing
end


"""
    haploimpute!(X, H, M, N, happair, hapscore, maxiters=1, tolfun=1e-3)

Haplotying of genotype matrix `X` from the pool of haplotypes `H` and impute
missing genotypes in `X` according to haplotypes.

# Input
* `X`: `n x p` nullable matrix. Each row is genotypes of an individual.
* `H`: `d x p` haplotype matrix. Each row is a haplotype.
* `M`: overwritten by `M[i, j] = 2dot(H[i, :], H[j, :]) + sumabs2(H[i, :]) +
    sumabs2(H[j, :])`.
* `N`: overwritten by `n x d` matrix `2XH'`.
* `happair`: optimal haplotype pair. `X[k, :] ≈ H[happair[k, 1], :] + H[happair[k, 2], :]`
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
* `maxiters`: number of MM iterations. Defaultis 1.
* `tolfun`: convergence tolerance of MM iterations. Default is 1e-3.
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
    initmissing!(X)
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

"""
    haploimpute!(X, H, width=500, maxiters=5, tolfun=1e-3)

Haplotying of genotype matrix `X` from a pool of haplotypes `H` and impute
missing genotypes in `X` according to haplotypes.

# Input
* `X`: `n x p` nullable matrix. Each row is genotypes of an individual.
* `H`: `d x p` haplotype matrix. Each row is a haplotype.
* `width`: width of the sliding window.
* `maxiters`: number of MM iterations. Defaultis 1.
* `tolfun`: convergence tolerance of MM iterations. Default is 1e-3.
"""
function haploimpute!(
    X::NullableMatrix,
    H::AbstractMatrix,
    width::Int     = 500,
    maxiters::Int  = 1,
    tolfun::Number = 1e-3,
    verbose::Bool  = true
    )

    people, snps, haplotypes = size(X, 1), size(X, 2), size(H, 1)
    # allocate working arrays
    M        = zeros(eltype(H), haplotypes, haplotypes)
    N        = zeros(promote_type(eltype(H), eltype(X.values)), people, haplotypes)
    happair  = zeros(Int, people), zeros(Int, people)
    hapscore = zeros(eltype(N), people)

    # no need for sliding window
    if snps ≤ 3width
        haploimpute!(X, H, M, N, happair, hapscore, maxiters, tolfun)
        fillgeno!(X.values, H, happair)
        return nothing
    end

    # allocate working arrays
    Xwork = X[:, 1:3width] # NullableMatrix
    Xw1   = view(Xwork.values, :, 1:width)
    Xwb1  = view(Xwork.isnull, :, 1:width)
    Xw23  = view(Xwork.values, :, (width + 1):3width)
    Xwb23 = view(Xwork.isnull, :, (width + 1):3width)
    Hwork = view(H, :, 1:3width)

    # number of windows
    windows = floor(Int, snps / width)

    # phase and impute window 1
    if verbose; println("Imputing SNPs 1:$width"); end
    haploimpute!(Xwork, Hwork, M, N, happair, hapscore, maxiters, tolfun)
    fill!(Xwb1, false)

    # first  1/3: ((w - 2) * width + 1):((w - 1) * width)
    # middle 1/3: ((w - 1) * width + 1):(w * width)
    # last   1/3:       (w * width + 1):((w + 1) * width)
    for w in 2:(windows - 1)
        if verbose
            println("Imputing SNPs $((w - 1) * width + 1):$(w * width)")
        end
        # overwrite first 1/3 by phased haplotypes
        H1    = view(H,        :, ((w - 2) * width + 1):((w - 1) * width))
        X1    = view(X.values, :, ((w - 2) * width + 1):((w - 1) * width))
        fillgeno!(X1, H1, happair)
        copy!(Xw1, X1)
        # refresh second and third 1/3 to original data
        X23   = view(X.values, :, ((w - 1) * width + 1):((w + 1) * width))
        Xb23  = view(X.isnull, :, ((w - 1) * width + 1):((w + 1) * width))
        copy!(Xw23, X23)
        copy!(Xwb23, Xb23)
        # phase + impute
        Hwork = view(H, :, ((w - 2) * width + 1):((w + 1) * width))
        haploimpute!(Xwork, Hwork, M, N, happair, hapscore, maxiters, tolfun)
    end

    # last window
    if verbose
        println("Imputing SNPs $((windows - 1) * width + 1):$snps")
    end
    Xwork = X[:, ((windows - 2) * width + 1):snps]
    Hwork = view(H,        :, ((windows - 2) * width + 1):snps)
    H1    = view(H,        :, ((windows - 2) * width + 1):((windows - 1) * width))
    X1    = view(X.values, :, ((windows - 2) * width + 1):((windows - 1) * width))
    fillgeno!(X1, H1, happair)
    copy!(Xw1, X1)
    H23   = view(H,            :, ((windows - 1) * width + 1):snps)
    X23   = view(X.values,     :, ((windows - 1) * width + 1):snps)
    Xw23  = view(Xwork.values, :, (width + 1):size(Xwork, 2))
    Xb23  = view(X.isnull,     :, ((windows - 1) * width + 1):snps)
    copy!(Xw23, X23)
    copy!(Xwb23, Xb23)
    haploimpute!(Xwork, Hwork, M, N, happair, hapscore, maxiters, tolfun)
    fillgeno!(X23, H23, happair)

    return nothing

end

"""
The version not using information from flanking windows.
"""
function haploimpute2!(
    X::NullableMatrix,
    H::AbstractMatrix,
    width::Int  = 500,
    maxiters::Int  = 1,
    tolfun::Number = 1e-3,
    verbose::Bool  = true
    )

    people, snps, haplotypes = size(X, 1), size(X, 2), size(H, 1)
    # allocate working arrays
    M        = zeros(eltype(H), haplotypes, haplotypes)
    N        = zeros(promote_type(eltype(H), eltype(X.values)), people, haplotypes)
    happair  = zeros(Int, people), zeros(Int, people)
    hapscore = zeros(eltype(N), people)

    # no need for sliding window
    if snps ≤ width
        haploimpute!(X, H, M, N, happair, hapscore, maxiters, tolfun)
        fillgeno!(X.values, H, happair)
        return nothing
    end

    # sliding window
    windows = floor(Int, snps / width)
    Xwork   = X[:, 1:width]
    for w in 1:(windows - 1)
        if verbose
            println("Imputing SNPs $((w - 1) * width + 1):$(w * width)")
        end
        Hwork = view(H,        :, ((w - 1) * width + 1):(w * width))
        Xdata = view(X.values, :, ((w - 1) * width + 1):(w * width))
        Xbool = view(X.isnull, :, ((w - 1) * width + 1):(w * width))
        #Xwork = NullableMatrix(Xdata, Xbool)
        #Xwork = view(X, :, ((w - 1) * width + 1):(w * width))
        copy!(Xwork.values, Xdata)
        copy!(Xwork.isnull, Xbool)
        haploimpute!(Xwork, Hwork, M, N, happair, hapscore, maxiters, tolfun)
        fillgeno!(Xdata, Hwork, happair)
    end

    # last window
    if verbose
        println("Imputing SNPs $((windows - 1) * width + 1):$snps")
    end
    Hwork = view(H, :, ((windows - 1) * width + 1):snps)
    Xwork = X[:, ((windows - 1) * width + 1):snps]
    Xdata = view(X.values, :, ((windows - 1) * width + 1):snps)
    haploimpute!(Xwork, Hwork, M, N, happair, hapscore, maxiters, tolfun)
    fillgeno!(Xdata, Hwork, happair)

    return nothing

end
