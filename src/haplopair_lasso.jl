######## THIS FILE IS THE SAME AS `haplotype_pair.jl`
######## except it solves the least squares objective using lasso penalty
######## || x - βH ||^2 + λ|| β ||_1 where we search λ so that || β ||_0 = 2 

function haplopair_lasso(
    X::AbstractMatrix,
    H::AbstractMatrix;
    keep::Int = 100
    )
    
    Xwork = zeros(Float32, size(X, 1), size(X, 2))
    Hwork = convert(Matrix{Float32}, H)
    initXfloat!(X, Xwork)

    p, n     = size(X)
    happairs = ones(Int, n), ones(Int, n)
    hapscore = zeros(Float32, n) # proportion of deviance explained for chosen model

    haplopair_lasso!(Xwork, Hwork, happairs, hapscore, keep)
    t1 = t2 = t3 = t4 = 0 # no haplotype rescreening

    return happairs, hapscore, t1, t2, t3, t4
end

"""
Finds λ so that || β ||_0 = 2 using a bijective search.
"""
function haplopair_lasso!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    happairs::Tuple{AbstractVector, AbstractVector},
    hapscore::AbstractVector
    )

    p, n, d = size(X, 1), size(X, 2), size(H, 2)

    for i in 1:n
        x = @view(X[:, i])
        
        # find λ so that only 2 β is non-zero
        λ_prev = 1.0
        λ_curr = 0.5
        while true
            path = fit(LassoPath, H, x, Normal(), IdentityLink(), λ=[λ_curr])

            nz = count(!iszero, path.coefs)
            if nz == 2
                # record answer
                hapscore[i] = path.pct_dev[1]
                happairs[1][i] = path.coefs.rowval[1]
                happairs[2][i] = path.coefs.rowval[2]
                break
            elseif nz > 2
                
            elseif nz < 2
                
            end
        end
    end

    return nothing
end
