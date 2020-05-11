
"""
    continue_haplotype(X, H, happair_prev, happair_next)

Find the optimal concatenated haplotypes from unordered haplotype pairs in two
consecutive windows.

# Input
* `X`: an `n` vector of genotypes with {0, 1, 2} entries
* `H`: an `n x d` reference panel of haplotypes with {0, 1} entries
* `happair_prev`: unordered haplotypes `(i, j)` in the first window
* `happair_next`: unordered haplotypes `(k, l)` in the second window

# Output
* `happair_next_optimal`: optimal ordered haplotypes in the second window
* `breakpt`: break points in the ordered haplotypes
"""
function continue_haplotype(
    X::AbstractVector,
    H::AbstractMatrix,
    happair_prev::Tuple{Int, Int},
    happair_next::Tuple{Int, Int}
    )

    i, j = happair_prev
    k, l = happair_next

    # both strands match
    if i == k && j == l
        return (k, l), (-1, -1)
    end

    if i == l && j == k
        return (l, k), (-1, -1)
    end

    # only one strand matches
    if i == k && j ≠ l
        breakpt, errors = search_breakpoint(X, H, i, (j, l))
        return (k, l), (-1, breakpt)
    elseif i == l && j ≠ k
        breakpt, errors = search_breakpoint(X, H, i, (j, k))
        return (l, k), (-1, breakpt)
    elseif j == k && i ≠ l
        breakpt, errors = search_breakpoint(X, H, j, (i, l))
        return (l, k), (breakpt, -1)
    elseif j == l && i ≠ k
        breakpt, errors = search_breakpoint(X, H, j, (i, k))
        return (k, l), (breakpt, -1)
    end

    # both strand mismatch
    breakpt1, errors1 = search_breakpoint(X, H, (i, k), (j, l))
    breakpt2, errors2 = search_breakpoint(X, H, (i, l), (j, k))
    if errors1 < errors2
        return (k, l), breakpt1
    else
        return (l, k), breakpt2
    end

    # width = round(Int, length(X) / 2) # must use round since last window width might not be integer
    # return (k, l), (width, width)
end

"""
    search_breakpoint(X, H, s1, s2)

Find the optimal break point between s2[1] and s2[2] in configuration
s1 | s2[1]
s1 | s2[2]
"""
function search_breakpoint(
    X::AbstractVector,
    H::AbstractMatrix,
    s1::Int,
    s2::Tuple{Int, Int}
    )

    n = length(X)
    # count number of errors if second haplotype is all from H[:, s2[2]]
    errors = 0
    for pos in 1:n
        if !ismissing(X[pos])
            errors += X[pos] ≠ H[pos, s1] + H[pos, s2[2]]
        end
    end
    bkpt_optim, err_optim = 0, errors

    # quick return if perfect match
    err_optim == 0 && return 0, 0

    # extend haplotype H[:, s2[1]] position by position
    @inbounds for bkpt in 1:n
        if !ismissing(X[bkpt]) && H[bkpt, s2[1]] ≠ H[bkpt, s2[2]]
            errors -= abs2(X[bkpt] - H[bkpt, s1] - H[bkpt, s2[2]])
            errors += abs2(X[bkpt] - H[bkpt, s1] - H[bkpt, s2[1]])
            if errors :: Int < err_optim
                bkpt_optim, err_optim = bkpt, errors
                # quick return if perfect match
                err_optim == 0 && return bkpt_optim, err_optim :: Int
            end
        end
    end

    return bkpt_optim, err_optim :: Int
end

"""
    search_breakpoint(X, H, s1, s2)

Find the optimal break point between s2[1] and s2[2] in configuration
s1[1] | s2[1]
s1[2] | s2[2]
"""
function search_breakpoint(
    X::AbstractVector,
    H::AbstractMatrix,
    s1::Tuple{Int, Int},
    s2::Tuple{Int, Int}
    )

    err_optim   = typemax(Int)
    bkpts_optim = (0, 0)

    # search over all combintations of break points in two strands
    @inbounds for bkpt1 in 0:length(X)

        # count number of errors if second haplotype is all from H[:, s2[2]]
        errors = 0
        for pos in 1:bkpt1
            if !ismissing(X[pos])
                errors += abs2(X[pos] - H[pos, s1[1]] - H[pos, s2[2]])
            end
        end
        for pos in (bkpt1 + 1):length(X)
            if !ismissing(X[pos])
                errors += abs2(X[pos] - H[pos, s1[2]] - H[pos, s2[2]])
            end
        end
        if errors :: Int < err_optim
            err_optim = errors
            bkpts_optim = (bkpt1, 0)

            # quick return if perfect match
            err_optim == 0 && return bkpts_optim, err_optim :: Int
        end

        # extend haplotype H[:, s2[1]] position by position
        for bkpt2 in 1:bkpt1
            if !ismissing(X[bkpt2]) && H[bkpt2, s2[2]] != H[bkpt2, s2[1]]
                errors -= abs2(X[bkpt2] - H[bkpt2, s1[1]] - H[bkpt2, s2[2]])
                errors += abs2(X[bkpt2] - H[bkpt2, s1[1]] - H[bkpt2, s2[1]])
                if errors :: Int < err_optim
                    err_optim = errors
                    bkpts_optim = (bkpt1, bkpt2)
                end
            end
        end
        for bkpt2 in (bkpt1 + 1):length(X)
            if !ismissing(X[bkpt2]) && H[bkpt2, s2[2]] != H[bkpt2, s2[1]]
                errors -= abs2(X[bkpt2] - H[bkpt2, s1[2]] - H[bkpt2, s2[2]])
                errors += abs2(X[bkpt2] - H[bkpt2, s1[2]] - H[bkpt2, s2[1]])
                if errors :: Int < err_optim
                    err_optim = errors
                    bkpts_optim = (bkpt1, bkpt2)
                    # quick return if perfect match
                    err_optim == 0 && return bkpts_optim, err_optim :: Int
                end
            end
        end
    end

    return bkpts_optim, err_optim :: Int
end

"""
    search_breakpoint_dp(X, H, s1, s2)

Find the optimal break point between s2[1] and s2[2] in configuration:
s1[1] | s2[1]
s1[2] | s2[2]

# Example:

    |----a-----|--c----|
    |--b----|----d-----|

If s1 = (a, b) and s2 = (c, d), then SNP in position 3 is (ab), position 10 is (ad), position 12 is (cd) and there are no (cb). 
Thus, snp is either (ab), (ad), (cb), or (cd). This order is assumed in vectors `sol_path`, `next_pair`, and `subtree_err`. 
"""
function search_breakpoint_dp(
    X::AbstractVector,
    H::AbstractMatrix,
    s1::Tuple{Int, Int},
    s2::Tuple{Int, Int}
    )

    # TODO: THIS FUNCTION IS NOT COMPLETE
    
    width = length(X)
    sol_path = zeros(Int, width)
    next_pair = [zeros(Int, 4) for i in 1:width]
    subtree_err = [zeros(Float64, 4) for i in 1:width]

    # iterate from 2nd to last snp, since last snp induces no error and connects to nothing
    for snp in Iterators.reverse(1:(width - 1))
        tree_err = [Inf for i in 1:4]
        next = zeros(Int, 4)

        for err in subtree_err[snp + 1]
            ab_err = abs2(X[snp] - H[snp, s1[1]] - H[snp, s2[1]]) + err
            ad_err = abs2(X[snp] - H[snp, s1[1]] - H[snp, s2[2]]) + err
            cb_err = abs2(X[snp] - H[snp, s1[2]] - H[snp, s2[1]]) + err
            cd_err = abs2(X[snp] - H[snp, s1[2]] - H[snp, s2[2]]) + err

            if ab_err < tree_err[1]
                tree_err[1] = ab_err
                next[1] = 1
            end
            if ad_err < tree_err[2]
                tree_err[2] = ad_err
                next[2] = 2
            end
            if cb_err < tree_err[3]
                tree_err[3] = cb_err
                next[3] = 3
            end
            if cd_err < tree_err[4]
                tree_err[4] = cd_err
                next[4] = 4
            end
        end

        # store current snp's best result for each of (ab), (ad), (cb), (cd)
        for i in 1:4
            subtree_err[i] = tree_err[i]
            next_pair[i] = next[i]
        end
    end

    return sol_path
end


function pair_error(Xi, h1, h2, cross)
    if cross > 1 
        return Inf
    else
        
    end
end


