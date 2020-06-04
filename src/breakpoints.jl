function continue_haplotype(
    X::AbstractVector,
    compressed_Hunique::CompressedHaplotypes,
    window::Int,
    happair_prev::Tuple{Int, Int},
    happair_next::Tuple{Int, Int}
    )

    # indices for complete reference panel
    i, j = happair_prev
    k, l = happair_next

    # both strands match
    if i == k && j == l
        return (k, l), (-1, -1)
    end
    if i == l && j == k
        return (l, k), (-1, -1)
    end

    Hprev = compressed_Hunique[window - 1].uniqueH
    Hcurr = compressed_Hunique[window].uniqueH

    # only one strand matches
    if i == k && j ≠ l
        iu = complete_idx_to_unique_idx(i, window - 1, compressed_Hunique)
        ku = complete_idx_to_unique_idx(k, window, compressed_Hunique)
        ju1 = complete_idx_to_unique_idx(j, window - 1, compressed_Hunique)
        ju2 = complete_idx_to_unique_idx(j, window, compressed_Hunique)
        lu1 = complete_idx_to_unique_idx(l, window - 1, compressed_Hunique)
        lu2 = complete_idx_to_unique_idx(l, window, compressed_Hunique)

        s1  = vcat(Hprev[:, iu], Hcurr[:, ku])
        s21 = vcat(Hprev[:, ju1], Hcurr[:, ju2])
        s22 = vcat(Hprev[:, lu1], Hcurr[:, lu2])

        breakpt, errors = search_breakpoint(X, s1, s21, s22)
        return (k, l), (-1, breakpt)
    elseif i == l && j ≠ k
        iu = complete_idx_to_unique_idx(i, window - 1, compressed_Hunique)
        lu = complete_idx_to_unique_idx(l, window, compressed_Hunique)
        ju1 = complete_idx_to_unique_idx(j, window - 1, compressed_Hunique)
        ju2 = complete_idx_to_unique_idx(j, window, compressed_Hunique)
        ku1 = complete_idx_to_unique_idx(k, window - 1, compressed_Hunique)
        ku2 = complete_idx_to_unique_idx(k, window, compressed_Hunique)

        s1  = vcat(Hprev[:, iu], Hcurr[:, lu])
        s21 = vcat(Hprev[:, ju1], Hcurr[:, ju2])
        s22 = vcat(Hprev[:, ku1], Hcurr[:, ku2])

        breakpt, errors = search_breakpoint(X, s1, s21, s22)
        return (l, k), (-1, breakpt)
    elseif j == k && i ≠ l
        ju = complete_idx_to_unique_idx(j, window - 1, compressed_Hunique)
        ku = complete_idx_to_unique_idx(k, window, compressed_Hunique)
        iu1 = complete_idx_to_unique_idx(i, window - 1, compressed_Hunique)
        iu2 = complete_idx_to_unique_idx(i, window, compressed_Hunique)
        lu1 = complete_idx_to_unique_idx(l, window - 1, compressed_Hunique)
        lu2 = complete_idx_to_unique_idx(l, window, compressed_Hunique)

        s1  = vcat(Hprev[:, ju], Hcurr[:, ku])
        s21 = vcat(Hprev[:, iu1], Hcurr[:, iu2])
        s22 = vcat(Hprev[:, lu1], Hcurr[:, lu2])

        breakpt, errors = search_breakpoint(X, s1, s21, s22)
        return (l, k), (breakpt, -1)
    elseif j == l && i ≠ k
        ju = complete_idx_to_unique_idx(j, window - 1, compressed_Hunique)
        lu = complete_idx_to_unique_idx(l, window, compressed_Hunique)
        iu1 = complete_idx_to_unique_idx(i, window - 1, compressed_Hunique)
        iu2 = complete_idx_to_unique_idx(i, window, compressed_Hunique)
        ku1 = complete_idx_to_unique_idx(k, window - 1, compressed_Hunique)
        ku2 = complete_idx_to_unique_idx(k, window, compressed_Hunique)

        s1  = vcat(Hprev[:, ju], Hcurr[:, lu])
        s21 = vcat(Hprev[:, iu1], Hcurr[:, iu2])
        s22 = vcat(Hprev[:, ku1], Hcurr[:, ku2])

        breakpt, errors = search_breakpoint(X, s1, s21, s22)
        return (k, l), (breakpt, -1)
    end

    # both strand mismatch
    # breakpt1, errors1 = search_breakpoint(X, H, (i, k), (j, l))
    # breakpt2, errors2 = search_breakpoint(X, H, (i, l), (j, k))
    # if errors1 < errors2
    #     return (k, l), breakpt1
    # else
    #     return (l, k), breakpt2
    # end
    return (k, l), (-1, -1)
end

"""
    search_single_breakpoint(X, s1, s21, s22)
Find the optimal break point between s21 and s22 in configuration
s1 | s21
s1 | s22
"""
function search_breakpoint(
    X::AbstractVector,
    s1::AbstractVector,
    s21::AbstractVector,
    s22::AbstractVector,
    )

    n = length(X)

    # count number of errors if second haplotype is all from s22
    errors = 0
    for pos in 1:n
        if !ismissing(X[pos])
            errors += X[pos] ≠ s1[pos] + s22[pos]
        end
    end
    bkpt_optim, err_optim = 0, errors

    # quick return if perfect match
    err_optim == 0 && return 0, 0

    # extend haplotype s21 position by position
    @inbounds for bkpt in 1:n
        if !ismissing(X[bkpt]) && s21[bkpt] ≠ s22[bkpt]
            errors -= X[bkpt] ≠ s1[bkpt] + s22[bkpt]
            errors += X[bkpt] ≠ s1[bkpt] + s21[bkpt]
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
            errors -= X[bkpt] ≠ H[bkpt, s1] + H[bkpt, s2[2]]
            errors += X[bkpt] ≠ H[bkpt, s1] + H[bkpt, s2[1]]
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
                errors += X[pos] ≠ H[pos, s1[1]] + H[pos, s2[2]]
            end
        end
        for pos in (bkpt1 + 1):length(X)
            if !ismissing(X[pos])
                errors += X[pos] ≠ H[pos, s1[2]] + H[pos, s2[2]]
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
            if !ismissing(X[bkpt2]) && H[bkpt2, s2[1]] ≠ H[bkpt2, s2[2]]
                errors -= X[bkpt2] ≠ H[bkpt2, s1[1]] + H[bkpt2, s2[2]]
                errors += X[bkpt2] ≠ H[bkpt2, s1[1]] + H[bkpt2, s2[1]]
                if errors :: Int < err_optim
                    err_optim = errors
                    bkpts_optim = (bkpt1, bkpt2)
                end
            end
        end
        for bkpt2 in (bkpt1 + 1):length(X)
            if !ismissing(X[bkpt2]) && H[bkpt2, s2[1]] ≠ H[bkpt2, s2[2]]
                errors -= X[bkpt2] ≠ H[bkpt2, s1[2]] + H[bkpt2, s2[2]]
                errors += X[bkpt2] ≠ H[bkpt2, s1[2]] + H[bkpt2, s2[1]]
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
