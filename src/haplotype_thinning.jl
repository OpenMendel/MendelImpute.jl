######## THIS FILE IS THE SAME AS `hpalotype_pair.jl`
######## except it solves the least squares objective on only "thining_factor" unique haplotypes.
######## where the top haplotypes are selected by minimizing dist(x, h). 

function haplopair_thin(
    X::AbstractMatrix,
    H::AbstractMatrix;
    keep::Int = 100
    )

    Xwork = zeros(Float32, size(X, 1), size(X, 2))
    Hwork = convert(Matrix{Float32}, H)
    initXfloat!(X, Xwork)

    p, n     = size(X)
    d        = size(H, 2)
    M        = zeros(Float32, d, d)
    N        = zeros(Float32, n, d)
    happairs = ones(Int, n), ones(Int, n)
    hapscore = zeros(Float32, n)

    t1, t2, t3 = haplopair!(Xwork, Hwork, M, N, happairs, hapscore, keep)
    t4 = 0 # no rescreening or 

    return happairs, hapscore, t1, t2, t3, t4
end

function haplopair!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    M::AbstractMatrix,
    N::AbstractMatrix,
    happairs::Tuple{AbstractVector, AbstractVector},
    hapscore::AbstractVector,
    keep::Int
    )

    p, n, d = size(X, 1), size(X, 2), size(H, 2)

    # compute distances between each column of H and each column of X
    t1 = @elapsed R = pairwise(Euclidean(), H, X, dims=2) # Rij = d(H[:, i], X[:, j])

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
    t3 = @elapsed haplopair!(happairs[1], happairs[2], hapscore, M, N, R, keep)

    # supplement the constant terms in objective
    t3 += @elapsed begin @inbounds for j in 1:n
            @simd for i in 1:p
                hapscore[j] += abs2(X[i, j])
            end
        end
    end

    return t1, t2, t3
end

function haplopair!(
    happair1::AbstractVector{Int},
    happair2::AbstractVector{Int},
    hapmin::AbstractVector{T},
    M::AbstractMatrix{T},
    N::AbstractMatrix{T},
    R::AbstractMatrix{T},
    keep::Int
    ) where T <: Real

    n, d = size(N)
    perm = zeros(Int, d)
    fill!(hapmin, typemax(T))
    keep > d && (keep = d)

    for i in 1:n
        # find top matching haplotypes for sample i
        partialsortperm!(perm, view(R, :, i), keep) # perm[1:keep] = hap indices that best matches xi 

        @inbounds for (idx1, k) in enumerate(perm)
            idx1 > keep && break
            for (idx2, j) in enumerate(perm)
                idx2 > keep && break
                k > j && continue # since M upper triangular
                
                score = M[k, j] - N[i, k] - N[i, j]
                if score < hapmin[i]
                    hapmin[i], happair1[i], happair2[i] = score, k, j
                end
            end
        end
    end

    return nothing
end
