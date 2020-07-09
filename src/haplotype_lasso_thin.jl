######## THIS FILE IS THE SAME AS `haplotype_thinning.jl`
######## except the top r haplotypes are selected via lasso

function haplopair_lasso_thin(
    X::AbstractMatrix,
    H::AbstractMatrix;
    r::Int = 100
    )

    p, n  = size(X)
    d     = size(H, 2)
    r > d && (r = d) #safety check

    Xwork = zeros(Float32, p, n)
    Hwork = convert(Matrix{Float32}, H)
    initXfloat!(X, Xwork)

    M   = zeros(Float32, r, r)
    idx = zeros(Int, d)
    Xi  = zeros(Float32, p)
    Ni  = zeros(Float32, r)
    Hk  = zeros(Float32, p, r)
    happairs = ones(Int, n), ones(Int, n)
    hapscore = zeros(Float32, n)

    # assemble (complete) N
    t2 = @elapsed begin
        Nt = zeros(Float32, d, n)
        mul!(Nt, Transpose(Hwork), Xwork)
        @simd for I in eachindex(Nt)
            Nt[I] *= 2
        end
    end

    # search all pairwise combination of top r haplotypes
    t1 = t3 = 0
    @inbounds for i in 1:n
        # find the top r haplotypes
        t1 += @elapsed begin
            @views sortperm!(idx, Nt[:, i]; alg = PartialQuickSort(r), rev = true)
        end

        # sync Hk, Xi, M, Ni
        t2 += @elapsed begin
            for k in 1:r
                col = idx[k]
                for j in 1:p
                    Hk[j, k] = Hwork[j, col]
                end
            end
            for j in eachindex(Xi)
                Xi[j] = Xwork[j, i]
            end
            update_M!(M, Hk)
            update_N!(Ni, Xi, Hk)
        end

        # search routine
        t3 += @elapsed begin
            hapscore[i], h1, h2 = haplopair!(M, Ni)
            happairs[1][i], happairs[2][i] = idx[h1], idx[h2]
        end

        # supply constant term in objective
        t2 += @elapsed hapscore[i] += dot(Xi, Xi)
    end

    t4 = 0 # no haplotype rescreening

    return happairs, hapscore, t1, t2, t3, t4
end 
