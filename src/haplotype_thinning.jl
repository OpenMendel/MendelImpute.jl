function haplopair_thin(
    X::AbstractMatrix,
    H::AbstractMatrix;
    keep::Int = 10
    )

    p, n = size(X)
    d    = size(H, 2)
    M    = zeros(Float32, keep, keep)
    N    = zeros(Float32, keep)
    perm = zeros(Int, d)
    Hk   = zeros(Float32, p, keep)
    Xi   = zeros(Float32, p)

    happairs = ones(Int, n), ones(Int, n)
    hapscore = zeros(Float32, n)

    # compute distances between each column of H and each column of X
    R = pairwise(Hamming(), H, X, dims=2) # Rij = Hamming(H[:, i], X[:, j])

    for i in 1:n
        # find top matching haplotypes for sample i
        partialsortperm!(perm, view(R, :, i), keep) # perm[1:keep] = hap indices that best matches xi 
        copyto!(Hk, @view(H[:, perm[1:keep]]))
        copyto!(Xi, @view(X[:, i]))
        update_M!(M, Hk)
        update_N!(N, Xi, Hk)

        # find best haplotype pair among the candidate pool
        hapscore[i], h1, h2 = haplopair!(M, N)
        happairs[1][i], happairs[2][i] = perm[h1], perm[h2]

        # supply constant term in objective
        hapscore[i] += dot(Xi, Xi)
    end

    return happairs, hapscore
end

function update_M!(M::AbstractMatrix, H::AbstractMatrix)
    d = size(M, 1)
    mul!(M, Transpose(H), H)
    for j in 1:d, i in 1:(j - 1) # off-diagonal
        M[i, j] = 2M[i, j] + M[i, i] + M[j, j]
    end
    for j in 1:d # diagonal
        M[j, j] *= 4
    end
    return nothing
end

function update_N!(N::AbstractVector, x::AbstractVector, H::AbstractMatrix)
    mul!(N, Transpose(H), x)
    @simd for I in eachindex(N)
        N[I] *= 2
    end
    return nothing
end

function haplopair!(
    M::AbstractMatrix{T},
    N::AbstractVector{T},
    ) where T <: Real

    d = length(N)
    hapmin, hap1, hap2 = typemax(T), 1, 1

    @inbounds for k in 1:d, j in 1:k
        score = M[j, k] - N[j] - N[k]
        if score < hapmin
            hapmin, hap1, hap2 = score, j, k
        end
    end

    return hapmin, hap1, hap2
end

n = 100
p = 100
d = 1000

H = convert(Matrix{Float32}, rand(0:1, p, d))
X = convert(Matrix{Float32}, rand(0:2, p, n))
keep = 10
perm = zeros(Int, d)
R = pairwise(Hamming(), H, X, dims=2) # Rij = Hamming(H[:, i], X[:, j])

i = 1
Ri = view(R, :, i)
partialsortperm!(perm, Ri, keep) # perm[1] = index of smallest Ri value




