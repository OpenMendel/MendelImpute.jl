######## THIS FILE IS THE SAME AS `haplotype_pair.jl`
######## except it solves the least squares objective on only `keep (default 100)` unique 
######## haplotypes, where the top haplotypes are selected by minimizing dist(x, h). 

function haplopair_thin_BLAS2!(
    X::AbstractMatrix,
    H::AbstractMatrix;
    alt_allele_freq::Union{Nothing, AbstractVector{Float32}} = nothing, 
    keep::Int = 100,
    # preallocated vectors
    happair1::AbstractVector = ones(Int, size(X, 2)),        # length n 
    happair2::AbstractVector = ones(Int, size(X, 2)),        # length n
    hapscore::AbstractVector = zeros(Float32, size(X, 2)),   # length n
    maxindx :: Vector{<:Integer} = Vector{Int}(undef, keep), # length keep
    maxgrad :: Vector{Float32}  = zeros(Float32, keep),      # length keep
    Xi :: AbstractVector{Float32} = zeros(Float32, size(H, 1)), # length p
    N  :: AbstractVector{Float32} = zeros(Float32, keep),       # length keep
    # preallocated matrices
    Hk    :: AbstractMatrix{Float32} = zeros(Float32, size(H, 1), keep),       # p × keep 
    M     :: AbstractMatrix{Float32} = zeros(Float32, keep, keep),             # keep × keep
    Xwork :: AbstractMatrix{Float32} = zeros(Float32, size(X, 1), size(X, 2)), # p × n
    Hwork :: AbstractMatrix{Float32} = convert(Matrix{Float32}, H),            # p × d
    R     :: AbstractMatrix{Float32} = zeros(Float32, size(H, 2), size(X, 2)), # d × p
    )

    p, n = size(X)
    d    = size(H, 2)
    keep > d && (keep = d) # safety check

    # reallocate matrices for last window
    if size(Hk, 1) != size(H, 1)
        Hk = zeros(Float32, size(H, 1), keep)
        Xi = zeros(Float32, size(H, 1))
        Xwork = zeros(Float32, p, n)

        # TODO: 
        # Hwork =
        # R = 
    end

    # sync Xwork and Hwork
    initXfloat!(Xwork, X)

    # compute distances between each column of H and each column of X
    t1 = @elapsed begin
        if !isnothing(alt_allele_freq)
            map!(x -> x < 0.5 ? 1 - x : x, alt_allele_freq, alt_allele_freq) # scale by 1 - p
            pairwise!(R, WeightedSqEuclidean(alt_allele_freq), Hwork, Xwork, dims=2)
        else 
            pairwise!(R, SqEuclidean(), Hwork, Xwork, dims=2) # Rij = || H[:, i] - X[:, j] ||²
        end 
    end

    t2 = t3 = 0
    @inbounds for i in 1:n
        # find top matching haplotypes for sample i
        t1 += @elapsed findbotr!(R, i, maxgrad, maxindx)

        # sync Hk, Xi, M, N
        t2 += @elapsed begin
            for k in 1:keep
                col = maxindx[k]
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
            happair1[i], happair2[i] = maxindx[h1], maxindx[h2]
        end

        # supply constant term in objective
        t2 += @elapsed hapscore[i] += dot(Xi, Xi)
    end

    t4 = 0 # no haplotype rescreening

    return t1, t2, t3, t4
end

"""
    findbotr!(A, col, maxval, index)

Find the smallest `length(index)` elements of column `A[:, col]`. `val` is 
filled with the found smallest `r` elements in sorted order and `idx` is filled
with their corresponding indices. 
"""
@inline function findbotr!(
    A      :: AbstractMatrix{T}, 
    col    :: Integer,
    minval :: AbstractVector{T}, 
    index  :: AbstractVector{<:Integer}
    ) where T <: Real
    fill!(minval, typemax(T))
    @inbounds for row in 1:size(A, 1)
        a = A[row, col]
        k = searchsortedlast(minval, a)
        if k < length(index)
            popinsert_rev!(minval, k+1, a)
            popinsert_rev!(index, k+1, row)
        end
    end
    nothing
end

"""
    popinsert_rev!(v, k, vk)

Move elements in `v[k:end-1]` to `v[k+1:end]` and insert `vk` at position `k` of vector `v`.
"""
@inline function popinsert_rev!(v::AbstractVector, k::Integer, vk)
    @inbounds for i in Iterators.reverse(k+1:length(v))
        v[i] = v[i-1]
    end
    v[k] = vk
    v
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

    @inbounds for k in 1:d
        Nk = N[k]
        for j in 1:k
            score = M[j, k] - N[j] - Nk
            if score < hapmin
                hapmin, hap1, hap2 = score, j, k
            end
        end
    end

    return hapmin, hap1, hap2
end

"""
Same as `haplopair_thin_BLAS2` but internally it uses only BLAS 3 calls and does not recompute `M`
for every sample. But this implies the search routine in `haplopair!` is not cache aware.  
"""
function haplopair_thin_BLAS3!(
    X::AbstractMatrix,
    H::AbstractMatrix;
    alt_allele_freq::Union{Nothing, AbstractVector{Float32}} = nothing, 
    keep::Int = 100,
    # N::ElasticArray{Float32} = ElasticArray{Float32}(undef, size(X, 2), size(H, 2)), # n × d
    happair1::AbstractVector = ones(Int, size(X, 2)),     # length n 
    happair2::AbstractVector = ones(Int, size(X, 2)),     # length n
    hapscore::AbstractVector = zeros(Float32, size(X, 2)) # length n
    )

    Xwork = zeros(Float32, size(X, 1), size(X, 2))
    Hwork = convert(Matrix{Float32}, H)
    initXfloat!(Xwork, X)

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

function haplopair_thin_BLAS3(
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
    t3 += @elapsed begin
        @inbounds for j in 1:n
            @simd for i in 1:p
                hapscore[j] += abs2(X[i, j])
            end
        end
    end

    return t1, t2, t3
end

function haplopair_thin_BLAS3!(
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
