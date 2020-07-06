######## THIS FILE IS THE SAME AS `haplotype_pair.jl`
######## except it solves the least squares objective on only `keep (default 100)` unique 
######## haplotypes, where the top haplotypes are selected by minimizing dist(x, h). 

function haplopair_thin(
    X::AbstractMatrix,
    H::AbstractMatrix;
    alt_allele_freq::Union{Nothing, AbstractVector{Float32}} = nothing, 
    keep::Int = 100
    )

    p, n = size(X)
    d    = size(H, 2)
    keep > d && (keep = d) # safety check

    Xwork = zeros(Float32, p, n)
    Hwork = convert(Matrix{Float32}, H)
    initXfloat!(X, Xwork)

    M    = zeros(Float32, keep, keep)
    N    = zeros(Float32, keep)
    perm = zeros(Int, d)
    Hk   = zeros(Float32, p, keep)
    Xi   = zeros(Float32, p)
    happairs = ones(Int, n), ones(Int, n)
    hapscore = zeros(Float32, n)
    t4 = 0 # no haplotype rescreening

    # compute distances between each column of H and each column of X
    t1 = @elapsed begin
        if !isnothing(alt_allele_freq)
            map!(x -> 1 / (2 * x * (1 - x)), alt_allele_freq, alt_allele_freq) # scale by 1 / 2p(1-p)
            R = pairwise(WeightedEuclidean(alt_allele_freq), Hwork, Xwork, dims=2)
        else 
            R = pairwise(Euclidean(), Hwork, Xwork, dims=2) # Rij = d(H[:, i], X[:, j])
        end 
    end

    t2 = t3 = 0
    @inbounds for i in 1:n
        # find top matching haplotypes for sample i
        partialsortperm!(perm, view(R, :, i), keep) # perm[1:keep] = hap indices that best matches xi 
        
        # sync Hk, Xi, M, N
        t2 += @elapsed begin
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
        t3 += @elapsed begin
            hapscore[i], h1, h2 = haplopair!(M, N)
            happairs[1][i], happairs[2][i] = perm[h1], perm[h2]
        end

        # supply constant term in objective
        t2 += @elapsed hapscore[i] += dot(Xi, Xi)
    end

    return happairs, hapscore, t1, t2, t3, t4
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

"""
Same as `haplopair_thin` but internally it uses only BLAS 3 calls and does not recompute `M`
for every sample. But this implies the search routine in `haplopair!` is not cache aware.  
"""
function haplopair_thin2(
    X::AbstractMatrix,
    H::AbstractMatrix;
    alt_allele_freq::Union{Nothing, AbstractVector{Float32}} = nothing, 
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

    t1, t2, t3 = haplopair_thin2!(Xwork, Hwork, M, N, happairs, hapscore, alt_allele_freq, keep)
    t4 = 0 # no haplotype rescreening

    return happairs, hapscore, t1, t2, t3, t4
end

function haplopair_thin2!(
    X::AbstractMatrix{Float32},
    H::AbstractMatrix{Float32},
    M::AbstractMatrix{Float32},
    N::AbstractMatrix{Float32},
    happairs::Tuple{AbstractVector, AbstractVector},
    hapscore::AbstractVector,
    alt_allele_freq::Union{Nothing, AbstractVector{Float32}},
    keep::Int
    )

    p, n, d = size(X, 1), size(X, 2), size(H, 2)

    # compute distances between each column of H and each column of X
    t1 = @elapsed begin
        if !isnothing(alt_allele_freq)
            map!(x -> 1 / (2 * x * (1 - x)), alt_allele_freq, alt_allele_freq) # scale by 1 / 2p(1-p)
            R = pairwise(WeightedEuclidean(alt_allele_freq), H, X, dims=2)
        else 
            R = pairwise(Euclidean(), H, X, dims=2) # Rij = d(H[:, i], X[:, j])
        end 
    end

    t2 = @elapsed begin 
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
    end

    # computational routine
    t3 = @elapsed haplopair_thin2!(happairs[1], happairs[2], hapscore, M, N, R, keep)

    # supplement the constant terms in objective
    t3 += @elapsed begin @inbounds for j in 1:n
            @simd for i in 1:p
                hapscore[j] += abs2(X[i, j])
            end
        end
    end

    return t1, t2, t3
end

function haplopair_thin2!(
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

    @inbounds for i in 1:n
        # find top matching haplotypes for sample i
        partialsortperm!(perm, view(R, :, i), keep) # perm[1:keep] = hap indices that best matches xi 

        for iter1 in 1:keep
            k = perm[iter1]
            for iter2 in 1:keep
                j = perm[iter2]
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
