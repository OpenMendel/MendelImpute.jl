"""
    haplopair!(happair, hapscore, M, N)

# Input
* `happair`: `n * 2` storage for haplotype index for each individual.
* `hapmin`: minimum offered by the optimal haplotype pair.
* `M`: `d x d` matrix `2HH'`, where `H` is the haplotype matrix with haplotypes
    in rows. Only the upper triangular part is used.
* `N`: `n x d` matrix `2XH'`, where `X` is the genotype matrix with individuals
    in rows.
"""
function haplopair!(
    happair::AbstractMatrix,
    hapmin::Vector,
    M::AbstractMatrix,
    N::AbstractMatrix
    )

    n, d = size(N)
    fill!(hapmin, typemax(eltype(hapmin)))
    for j in 1:d, i in 1:j
        mij = M[i, j] + M[i, i] + M[j, j]
        # loop over individuals
        @inbounds @simd for k in 1:n
            score = mij - N[k, i] - N[k, j]
            if score < hapmin[k]
                hapmin[k]   = score
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
* `H`: `d x p` haplotype matrix. Each row is a haplotype. s

# Output
* `happair`: haplotyping result. `X[k, :] ≈ H[happair[k, 1], :] + H[happair[k, 2], :]`
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
"""
function haplopair(
    X::AbstractMatrix,
    H::AbstractMatrix
    )

    n, d = size(X, 1), size(H, 1)
    M = A_mul_Bt(H, H)
    M .*= 2
    N = A_mul_Bt(X, H)
    N .*= 2
    happair  = zeros(Int, n, 2)
    hapscore = zeros(eltype(N), n)
    haplopair!(happair, hapscore, M, N)
    for i in 1:n
        hapscore[i] += sumabs2(view(X, i, :))
    end

end
