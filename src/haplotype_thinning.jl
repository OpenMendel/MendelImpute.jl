function haplopair_thin(
    X::AbstractMatrix,
    H::AbstractMatrix;
    keep::Int = 100
    )

    Xwork = zeros(Float32, size(X, 1), size(X, 2))
    Hwork = convert(Matrix{Float32}, H)
    initXfloat!(X, Xwork)

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
    t1 = @elapsed R = pairwise(Euclidean(), Hwork, Xwork, dims=2) # Rij = d(H[:, i], X[:, j])
    # t1 = @elapsed R = pairwise(SqEuclidean(), Hwork, Xwork, dims=2) # Rij = d(H[:, i], X[:, j])
    # t1 = @elapsed R = pairwise(Hamming(), Hwork, Xwork, dims=2) # Rij = d(H[:, i], X[:, j])

    t1 = t2 = t3 = 0
    for i in 1:n
        # find top matching haplotypes for sample i
        partialsortperm!(perm, view(R, :, i), keep) # perm[1:keep] = hap indices that best matches xi 
        
        # sync Hk, Xi, M, N
        t1 += @elapsed begin
            for k in 1:keep
                col = perm[k]
                for j in 1:p
                    Hk[j, k] = Hwork[j, col]
                end
            end
            for j in eachindex(Xi)
                Xi[j] = Xwork[j, i]
            end
            update_M!(M, Hk)
            update_N!(N, Xi, Hk)
        end

        # computational routine
        t2 += @elapsed begin
            hapscore[i], h1, h2 = haplopair!(M, N)
            happairs[1][i], happairs[2][i] = perm[h1], perm[h2]
        end

        # supply constant term in objective
        t3 += @elapsed hapscore[i] += dot(Xi, Xi)
    end

    return happairs, hapscore, t1, t2, t3
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

# function haplopair_thin2(
#     X::AbstractMatrix,
#     H::AbstractMatrix;
#     keep::Int = 100
#     )

#     Xwork = zeros(Float32, size(X, 1), size(X, 2))
#     Hwork = convert(Matrix{Float32}, H)
#     initXfloat!(X, Xwork)

#     p, n     = size(X)
#     d        = size(H, 2)
#     M        = zeros(Float32, d, d)
#     N        = zeros(Float32, n, d)
#     happairs = ones(Int, n), ones(Int, n)
#     hapscore = zeros(Float32, n)

#     t1, t2, t3 = haplopair!(Xwork, Hwork, M, N, happairs, hapscore, keep)

#     return happairs, hapscore, t1, t2, t3
# end

# function haplopair!(
#     X::AbstractMatrix,
#     H::AbstractMatrix,
#     M::AbstractMatrix,
#     N::AbstractMatrix,
#     happairs::Tuple{AbstractVector, AbstractVector},
#     hapscore::AbstractVector,
#     keep::Int
#     )

#     p, n, d = size(X, 1), size(X, 2), size(H, 2)
#     perm = zeros(Int, d)

#     # compute distances between each column of H and each column of X
#     t1 = @elapsed R = pairwise(Euclidean(), H, X, dims=2) # Rij = d(H[:, i], X[:, j])

#     # assemble M (upper triangular only)
#     t1 += @elapsed begin 
#         mul!(M, Transpose(H), H)
#         for j in 1:d, i in 1:(j - 1) # off-diagonal
#             M[i, j] = 2M[i, j] + M[i, i] + M[j, j]
#         end
#         for j in 1:d # diagonal
#             M[j, j] *= 4
#         end

#         # assemble N
#         mul!(N, Transpose(X), H)
#         @simd for I in eachindex(N)
#             N[I] *= 2
#         end
#     end

#     # computational routine
#     t2 = @elapsed haplopair!(happairs[1], happairs[2], hapscore, M, N, R, perm, keep)

#     # supplement the constant terms in objective
#     t3 = @elapsed begin @inbounds for j in 1:n
#             @simd for i in 1:p
#                 hapscore[j] += abs2(X[i, j])
#             end
#         end
#     end

#     return t1, t2, t3
# end

# function haplopair!(
#     happair1::AbstractVector{Int},
#     happair2::AbstractVector{Int},
#     hapmin::AbstractVector{T},
#     M::AbstractMatrix{T},
#     N::AbstractMatrix{T},
#     R::AbstractMatrix{T},
#     perm::AbstractVector{Int},
#     keep::Int
#     ) where T <: Real

#     n, d = size(N)
#     fill!(hapmin, typemax(T))

#     for i in 1:n
#         # find top matching haplotypes for sample i
#         partialsortperm!(perm, view(R, :, i), keep) # perm[1:keep] = hap indices that best matches xi 
        
#         @inbounds for (idx1, k) in enumerate(perm)
#             idx1 > keep && break
#             for (idx2, j) in enumerate(perm)
#                 idx2 > keep && break
#                 j > k && continue

#                 score = M[k, j] - N[i, k] - N[i, j]
#                 if score < hapmin[i]
#                     hapmin[i], happair1[i], happair2[i] = score, k, j
#                 end
#             end
#         end
#     end    

#     return nothing
# end
