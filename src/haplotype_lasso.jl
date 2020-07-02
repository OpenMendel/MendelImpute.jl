######## THIS FILE IS THE SAME AS `haplotype_pair.jl`
######## except it solves the least squares objective using lasso penalty
######## || x - βH ||^2 + λ|| β ||_1 where we search λ so that || β ||_0 = 2 

function haplopair_lasso(
    X::AbstractMatrix,
    H::AbstractMatrix
    )
    
    Xwork = zeros(Float32, size(X, 1), size(X, 2))
    Hwork = convert(Matrix{Float32}, H)
    initXfloat!(X, Xwork)

    p, n     = size(X)
    happairs = ones(Int, n), ones(Int, n)
    hapscore = zeros(Float32, n) # proportion of deviance explained for chosen model

    haplopair_lasso!(Xwork, Hwork, happairs, hapscore)
    t1 = t2 = t3 = t4 = 0 # no haplotype rescreening

    return happairs, hapscore, t1, t2, t3, t4
end

"""
Finds λ so that only 2 β is non-zero using a bisection search.
"""
function haplopair_lasso!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    happairs::Tuple{AbstractVector, AbstractVector},
    hapscore::AbstractVector
    )

    p, n, d = size(X, 1), size(X, 2), size(H, 2)

    for i in 1:n
        x = X[:, i]
        
        # find λ so that only 2 β is non-zero
        λtop = 1.0
        λlow = 0.0
        linesearch = 0
        while true
            λmid = (λtop + λlow) / 2
            path = fit(LassoPath, H, x, Normal(), IdentityLink(), λ=[λmid])

            nz = count(!iszero, path.coefs)
            if nz == 2 || (linesearch > 10 && nz > 1) # dangerous hack to bypass cases where nz is never 2
                # record answer
                hapscore[i] = path.pct_dev[1]
                happairs[1][i] = path.coefs.rowval[1]
                happairs[2][i] = path.coefs.rowval[2]
                # println("λ = $λmid")
                break
            else
                # bisection search
                nz < 2 ? (λtop = λmid) : (λlow = λmid)
                linesearch += 1
            end
            # println("λtop = $λtop, λlow = $λlow, λmid = $λmid, nz = $nz")
        end
    end

    return nothing
end
