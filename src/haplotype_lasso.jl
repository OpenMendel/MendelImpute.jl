######## THIS FILE IS THE SAME AS `haplotype_pair.jl`
######## except it solves the least squares objective using a 
######## heuristic stepwise search

function haplopair_lasso(
    X::AbstractMatrix,
    H::AbstractMatrix;
    r::Int = 1
    )
    
    p, n  = size(X)
    d     = size(H, 2)
    r > d && (r = d) #safety check

    Xwork = zeros(Float32, p, n)
    Hwork = convert(Matrix{Float32}, H)
    initXfloat!(X, Xwork)

    happairs = ones(Int, n), ones(Int, n)
    hapscore = zeros(Float32, n)

    # assemble M (symmetric)
    stamp = time()
    M = zeros(Float32, d, d)
    mul!(M, Transpose(Hwork), Hwork)
    for j in 1:d, i in 1:(j - 1) # off-diagonal
        M[i, j] = 2M[i, j] + M[i, i] + M[j, j]
    end
    for j in 1:d # diagonal
        M[j, j] *= 4
    end
    LinearAlgebra.copytri!(M, 'U')
    t2 = time() - stamp

    # assemble N
    stamp = time()
    Nt = zeros(Float32, d, n)
    mul!(Nt, Transpose(Hwork), Xwork)
    @simd for I in eachindex(Nt)
        Nt[I] *= 2
    end
    t2 += time() - stamp

    # computational routine
    stamp = time()
    if r == 1
        haplopair_stepwise!(happairs[1], happairs[2], hapscore, M, Nt)
    else 
        haplopair_topr!(happairs[1], happairs[2], hapscore, M, Nt, r = r)
    end
    t3 = time() - stamp

    # supplement the constant terms in objective
    stamp = time()
    @inbounds for j in 1:n
        @simd for i in 1:p
            hapscore[j] += abs2(Xwork[i, j])
        end
    end
    t2 += time() - stamp

    t1 = t4 = 0 # no haplotype rescreening or computing dist(X, H)

    return happairs, hapscore, t1, t2, t3, t4
end 

function haplopair_stepwise!(
    happair1 :: AbstractVector{<:Integer},
    happair2 :: AbstractVector{<:Integer},
    hapmin   :: AbstractVector{T},
    M        :: AbstractMatrix{T}, # d x d
    Nt       :: AbstractMatrix{T}, # d x n
    ) where T <: Real
    d, n = size(Nt)
    fill!(hapmin, typemax(T))
    @inbounds for k in 1:n
        # find the first haplotype
        gmax, i1 = Nt[1, k], 1
        for i in 2:d
            g = Nt[i, k]
            if g > gmax
                gmax = g
                i1   = i
            end
        end
        happair1[k] = i1
        # find the optimal second haplotype given i1 using LS criterion
        for i in 1:d
            score = M[i, i1] - Nt[i, k]
            if score < hapmin[k]
                hapmin[k] = score
                happair2[k] = i
            end
        end
        hapmin[k] -= Nt[i1, k]
    end
    return nothing
end

# same as haplopair_stepwise! but searches top `r` haplotypes with largest gradient
function haplopair_topr!(
    happair1 :: AbstractVector{<:Integer},
    happair2 :: AbstractVector{<:Integer},
    hapmin   :: AbstractVector{T},
    M        :: AbstractMatrix{T}, # d x d
    Nt       :: AbstractMatrix{T}; # d x n
    r        :: Integer = 5
    ) where T <: Real
    d, n = size(Nt)
    fill!(hapmin, typemax(T))
    idx = Vector{Int}(undef, d)
    @inbounds for k in 1:n
        # find the top r haplotypes
        @views sortperm!(idx, Nt[:, k]; alg = PartialQuickSort(r), rev = true)
        # for each top haplotype, find the optimal second one
        for riter in 1:r
            i1 = idx[riter]
            Nt_i1k = Nt[i1, k]
            for i in 1:d
                score = M[i, i1] - Nt_i1k - Nt[i, k]
                if score < hapmin[k]
                    hapmin[k] = score
                    happair1[k] = i1
                    happair2[k] = i
                end
            end
        end
    end
    return nothing
end
