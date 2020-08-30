###### This file is part of the MendelImpute.jl package.
###### It contains relevant code to search breakpoints. 
###### Currently code for double breakpoint search is broken and commented out. 

"""
    continue_haplotype(X, compressed_Hunique, window, happair_prev, happair_next)

Searches the breakpoint between `happair_prev` and `happair_next`, if there is
one. Currently breakpoint is set to the middle if there is double breakpoints. 

# Arguments:
- `X`: Genotype vector spanning 2 windows. 
- `compressed_Hunique`: A `CompressedHaplotypes` keeping track of unique
    haplotypes for each window and some other information
- `window`: The current window being considered. It is the 2nd window of `X`. 
- `happair_prev`: Haplotype pair for first window
- `happair_next`: Haplotype pair for second window
- `phased`: Boolean indicating whether `happair_prev` and `happair_next` have
    been phased. If `false`, will try 2 different orientations. 

# Output
- `happair_next`: If `phase = false`, this is equal to `happair_next` for input. 
    Otherwise it may be flipped
- `bkpt`: Tuple of integer indicating optimal breakpoint starting from the 
    previous window. `-1` indicates no breakpoints and `-2` indicates double 
    breakpoint. 
"""
function continue_haplotype(
    X::AbstractVector,
    compressed_Hunique::CompressedHaplotypes,
    window::Int,
    happair_prev::Tuple,
    happair_next::Tuple;
    phased::Bool = false
    )

    # indices for complete reference panel
    i, j = happair_prev
    k, l = happair_next

    # both strands match
    if i == k && j == l
        return (k, l), (-1, -1)
    end
    if i == l && j == k && !phased
        return (l, k), (-1, -1)
    end

    # unique haplotypes with only typed snps
    Hprev = compressed_Hunique.CW_typed[window - 1].uniqueH
    Hcurr = compressed_Hunique.CW_typed[window].uniqueH

    # only one strand matches
    if i == k && j ≠ l
        iu = complete_idx_to_unique_typed_idx(i, window - 1, compressed_Hunique)
        ku = complete_idx_to_unique_typed_idx(k, window, compressed_Hunique)
        ju1 = complete_idx_to_unique_typed_idx(j, window - 1, compressed_Hunique)
        ju2 = complete_idx_to_unique_typed_idx(j, window, compressed_Hunique)
        lu1 = complete_idx_to_unique_typed_idx(l, window - 1, compressed_Hunique)
        lu2 = complete_idx_to_unique_typed_idx(l, window, compressed_Hunique)

        # TODO: use `LazyArrays.jl` for `vcat` 
        s1  = vcat(Hprev[:, iu], Hcurr[:, ku])
        s21 = vcat(Hprev[:, ju1], Hcurr[:, ju2])
        s22 = vcat(Hprev[:, lu1], Hcurr[:, lu2])

        breakpt, errors = search_breakpoint(X, s1, s21, s22)
        return (k, l), (-1, breakpt)
    elseif i == l && j ≠ k && !phased
        iu = complete_idx_to_unique_typed_idx(i, window - 1, compressed_Hunique)
        lu = complete_idx_to_unique_typed_idx(l, window, compressed_Hunique)
        ju1 = complete_idx_to_unique_typed_idx(j, window - 1, compressed_Hunique)
        ju2 = complete_idx_to_unique_typed_idx(j, window, compressed_Hunique)
        ku1 = complete_idx_to_unique_typed_idx(k, window - 1, compressed_Hunique)
        ku2 = complete_idx_to_unique_typed_idx(k, window, compressed_Hunique)

        # TODO: use `LazyArrays.jl` for `vcat` 
        s1  = vcat(Hprev[:, iu], Hcurr[:, lu])
        s21 = vcat(Hprev[:, ju1], Hcurr[:, ju2])
        s22 = vcat(Hprev[:, ku1], Hcurr[:, ku2])

        breakpt, errors = search_breakpoint(X, s1, s21, s22)
        return (l, k), (-1, breakpt)
    elseif j == k && i ≠ l && !phased
        ju = complete_idx_to_unique_typed_idx(j, window - 1, compressed_Hunique)
        ku = complete_idx_to_unique_typed_idx(k, window, compressed_Hunique)
        iu1 = complete_idx_to_unique_typed_idx(i, window - 1, compressed_Hunique)
        iu2 = complete_idx_to_unique_typed_idx(i, window, compressed_Hunique)
        lu1 = complete_idx_to_unique_typed_idx(l, window - 1, compressed_Hunique)
        lu2 = complete_idx_to_unique_typed_idx(l, window, compressed_Hunique)

        # TODO: use `LazyArrays.jl` for `vcat` 
        s1  = vcat(Hprev[:, ju], Hcurr[:, ku])
        s21 = vcat(Hprev[:, iu1], Hcurr[:, iu2])
        s22 = vcat(Hprev[:, lu1], Hcurr[:, lu2])

        breakpt, errors = search_breakpoint(X, s1, s21, s22)
        return (l, k), (breakpt, -1)
    elseif j == l && i ≠ k
        ju = complete_idx_to_unique_typed_idx(j, window - 1, compressed_Hunique)
        lu = complete_idx_to_unique_typed_idx(l, window, compressed_Hunique)
        iu1 = complete_idx_to_unique_typed_idx(i, window - 1, compressed_Hunique)
        iu2 = complete_idx_to_unique_typed_idx(i, window, compressed_Hunique)
        ku1 = complete_idx_to_unique_typed_idx(k, window - 1, compressed_Hunique)
        ku2 = complete_idx_to_unique_typed_idx(k, window, compressed_Hunique)

        # TODO: use `LazyArrays.jl` for `vcat` 
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
    return (k, l), (-2, -2)
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
    length(s1) == length(s21) == length(s22) == n ||error("search_breakpoint:" * 
        " all vectors should have same length but length X = $n, s1 =" * 
        " $(length(s1)), s21 = $(length(s21)), s22 = $(length(s22))")

    # count number of errors if second haplotype is all from s22
    errors = 0
    @inbounds for pos in 1:n
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



# """
#     search_breakpoint(X, H, s1, s2)
# Find the optimal break point between s2[1] and s2[2] in configuration
# s1[1] | s2[1]
# s1[2] | s2[2]
# """
# function search_breakpoint(
#     X::AbstractVector,
#     H::AbstractMatrix,
#     s1::Tuple{Int, Int},
#     s2::Tuple{Int, Int}
#     )

#     err_optim   = typemax(Int)
#     bkpts_optim = (0, 0)

#     # search over all combintations of break points in two strands
#     @inbounds for bkpt1 in 0:length(X)

#         # count number of errors if second haplotype is all from H[:, s2[2]]
#         errors = 0
#         for pos in 1:bkpt1
#             if !ismissing(X[pos])
#                 errors += X[pos] ≠ H[pos, s1[1]] + H[pos, s2[2]]
#             end
#         end
#         for pos in (bkpt1 + 1):length(X)
#             if !ismissing(X[pos])
#                 errors += X[pos] ≠ H[pos, s1[2]] + H[pos, s2[2]]
#             end
#         end
#         if errors :: Int < err_optim
#             err_optim = errors
#             bkpts_optim = (bkpt1, 0)

#             # quick return if perfect match
#             err_optim == 0 && return bkpts_optim, err_optim :: Int
#         end

#         # extend haplotype H[:, s2[1]] position by position
#         for bkpt2 in 1:bkpt1
#             if !ismissing(X[bkpt2]) && H[bkpt2, s2[1]] ≠ H[bkpt2, s2[2]]
#                 errors -= X[bkpt2] ≠ H[bkpt2, s1[1]] + H[bkpt2, s2[2]]
#                 errors += X[bkpt2] ≠ H[bkpt2, s1[1]] + H[bkpt2, s2[1]]
#                 if errors :: Int < err_optim
#                     err_optim = errors
#                     bkpts_optim = (bkpt1, bkpt2)
#                 end
#             end
#         end
#         for bkpt2 in (bkpt1 + 1):length(X)
#             if !ismissing(X[bkpt2]) && H[bkpt2, s2[1]] ≠ H[bkpt2, s2[2]]
#                 errors -= X[bkpt2] ≠ H[bkpt2, s1[2]] + H[bkpt2, s2[2]]
#                 errors += X[bkpt2] ≠ H[bkpt2, s1[2]] + H[bkpt2, s2[1]]
#                 if errors :: Int < err_optim
#                     err_optim = errors
#                     bkpts_optim = (bkpt1, bkpt2)
#                     # quick return if perfect match
#                     err_optim == 0 && return bkpts_optim, err_optim :: Int
#                 end
#             end
#         end
#     end

#     return bkpts_optim, err_optim :: Int
# end
