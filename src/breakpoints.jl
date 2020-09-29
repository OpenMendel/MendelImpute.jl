###### This file is part of the MendelImpute.jl package.
###### It contains relevant code to search breakpoints. 

"""
    continue_haplotype(X, compressed_Hunique, window, happair_prev, happair_next)

Searches the breakpoint between `happair_prev` and `happair_next`, if there is
one. 

# Arguments:
- `X`: Genotype vector spanning 2 windows. 
- `compressed_Hunique`: A `CompressedHaplotypes` keeping track of unique
    haplotypes for each window and some other information
- `window`: The current window being considered. It is the 2nd window of `X`. 
- `happair_prev`: Haplotype pair for first window
- `happair_next`: Haplotype pair for second window

# Optional arguments:
- `phased`: Boolean indicating whether `happair_prev` and `happair_next` have
    been phased. If `false`, will try 2 different orientations. 
- `search_double_bkpts`: If `true`, will search double breakpoints. If false 
    the output `bkpt = (-2, -2)`

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
    phased::Bool = false,
    search_double_bkpts::Bool = false
    )

    # indices for complete reference panel
    i, j = happair_prev
    k, l = happair_next

    # both strands match
    if i == k && j == l
        return (k, l), (-1, -1), 0
    end
    if i == l && j == k && !phased
        return (l, k), (-1, -1), 0
    end

    # get unique haplotypes with only typed snps
    overlap = compressed_Hunique.overlap
    if overlap
        # create views on non-overlapping regions
        rprev = nonoverlap_range(compressed_Hunique, window - 1)
        rcurr = nonoverlap_range(compressed_Hunique, window)
        Hprev = view(compressed_Hunique.CW_typed[window - 1].uniqueH, rprev, :)
        Hcurr = view(compressed_Hunique.CW_typed[window].uniqueH, rcurr, :)
    else
        Hprev = compressed_Hunique.CW_typed[window - 1].uniqueH
        Hcurr = compressed_Hunique.CW_typed[window].uniqueH
    end

    # only one strand matches
    if i == k && j ≠ l
        iu = complete_idx_to_unique_typed_idx(i, window - 1, compressed_Hunique)
        ku = complete_idx_to_unique_typed_idx(k, window, compressed_Hunique)
        ju1 = complete_idx_to_unique_typed_idx(j, window - 1, compressed_Hunique)
        ju2 = complete_idx_to_unique_typed_idx(j, window, compressed_Hunique)
        lu1 = complete_idx_to_unique_typed_idx(l, window - 1, compressed_Hunique)
        lu2 = complete_idx_to_unique_typed_idx(l, window, compressed_Hunique)

        # lazy concatenation 
        s1  = ApplyArray(vcat, view(Hprev, :, iu),  view(Hcurr, :, ku))
        s21 = ApplyArray(vcat, view(Hprev, :, ju1), view(Hcurr, :, ju2))
        s22 = ApplyArray(vcat, view(Hprev, :, lu1), view(Hcurr, :, lu2))

        breakpt, errors = search_breakpoint(X, s1, s21, s22)
        return (k, l), (-1, breakpt), errors
    elseif i == l && j ≠ k && !phased
        iu = complete_idx_to_unique_typed_idx(i, window - 1, compressed_Hunique)
        lu = complete_idx_to_unique_typed_idx(l, window, compressed_Hunique)
        ju1 = complete_idx_to_unique_typed_idx(j, window - 1, compressed_Hunique)
        ju2 = complete_idx_to_unique_typed_idx(j, window, compressed_Hunique)
        ku1 = complete_idx_to_unique_typed_idx(k, window - 1, compressed_Hunique)
        ku2 = complete_idx_to_unique_typed_idx(k, window, compressed_Hunique)

        # lazy concatenation 
        s1  = ApplyArray(vcat, view(Hprev, :, iu),  view(Hcurr, :, lu))
        s21 = ApplyArray(vcat, view(Hprev, :, ju1), view(Hcurr, :, ju2))
        s22 = ApplyArray(vcat, view(Hprev, :, ku1), view(Hcurr, :, ku2))

        breakpt, errors = search_breakpoint(X, s1, s21, s22)
        return (l, k), (-1, breakpt), errors
    elseif j == k && i ≠ l && !phased
        ju = complete_idx_to_unique_typed_idx(j, window - 1, compressed_Hunique)
        ku = complete_idx_to_unique_typed_idx(k, window, compressed_Hunique)
        iu1 = complete_idx_to_unique_typed_idx(i, window - 1, compressed_Hunique)
        iu2 = complete_idx_to_unique_typed_idx(i, window, compressed_Hunique)
        lu1 = complete_idx_to_unique_typed_idx(l, window - 1, compressed_Hunique)
        lu2 = complete_idx_to_unique_typed_idx(l, window, compressed_Hunique)
        
        # lazy concatenation 
        s1  = ApplyArray(vcat, view(Hprev, :, ju),  view(Hcurr, :, ku))
        s21 = ApplyArray(vcat, view(Hprev, :, iu1), view(Hcurr, :, iu2))
        s22 = ApplyArray(vcat, view(Hprev, :, lu1), view(Hcurr, :, lu2))

        breakpt, errors = search_breakpoint(X, s1, s21, s22)
        return (l, k), (breakpt, -1), errors
    elseif j == l && i ≠ k
        ju = complete_idx_to_unique_typed_idx(j, window - 1, compressed_Hunique)
        lu = complete_idx_to_unique_typed_idx(l, window, compressed_Hunique)
        iu1 = complete_idx_to_unique_typed_idx(i, window - 1, compressed_Hunique)
        iu2 = complete_idx_to_unique_typed_idx(i, window, compressed_Hunique)
        ku1 = complete_idx_to_unique_typed_idx(k, window - 1, compressed_Hunique)
        ku2 = complete_idx_to_unique_typed_idx(k, window, compressed_Hunique)

        # lazy concatenation 
        s1  = ApplyArray(vcat, view(Hprev, :, ju),  view(Hcurr, :, lu))
        s21 = ApplyArray(vcat, view(Hprev, :, iu1), view(Hcurr, :, iu2))
        s22 = ApplyArray(vcat, view(Hprev, :, ku1), view(Hcurr, :, ku2))

        breakpt, errors = search_breakpoint(X, s1, s21, s22)
        return (k, l), (breakpt, -1), errors
    end

    # both strand mismatch
    if search_double_bkpts
        iu1 = complete_idx_to_unique_typed_idx(i, window - 1, compressed_Hunique)
        iu2 = complete_idx_to_unique_typed_idx(i, window, compressed_Hunique)
        ku1 = complete_idx_to_unique_typed_idx(k, window - 1, compressed_Hunique)
        ku2 = complete_idx_to_unique_typed_idx(k, window, compressed_Hunique)
        ju1 = complete_idx_to_unique_typed_idx(j, window - 1, compressed_Hunique)
        ju2 = complete_idx_to_unique_typed_idx(j, window, compressed_Hunique)
        lu1 = complete_idx_to_unique_typed_idx(l, window - 1, compressed_Hunique)
        lu2 = complete_idx_to_unique_typed_idx(l, window, compressed_Hunique)
        s11 = ApplyArray(vcat, view(Hprev, :, iu1), view(Hcurr, :, iu2))
        s12 = ApplyArray(vcat, view(Hprev, :, ku1), view(Hcurr, :, ku2))
        s21 = ApplyArray(vcat, view(Hprev, :, ju1), view(Hcurr, :, ju2))
        s22 = ApplyArray(vcat, view(Hprev, :, lu1), view(Hcurr, :, lu2))

        breakpt1, errors1 = search_breakpoint(X, s11, s12, s21, s22)
        return (k, l), breakpt1, errors1

        # No need for 2nd double bkpt config, since window-by-window intersection
        # already took care of phasing  
        # breakpt1, errors1 = search_breakpoint(X, H, (i, k), (j, l))
        # breakpt2, errors2 = search_breakpoint(X, H, (i, l), (j, k))
        # if errors1 < errors2
        #     return (k, l), breakpt1, errors1
        # else
        #     return (l, k), breakpt2, errors2
        # end
    end

    return (k, l), (-2, -2), 0
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

"""
    search_breakpoint(X, s11, s12, s21, s22)

Find the optimal break point between s11 and s12 as well as between s21 and s22
in configuration
s11 | s21
s12 | s22
"""
function search_breakpoint(
    X::AbstractVector,
    s11::AbstractVector,
    s12::AbstractVector,
    s21::AbstractVector,
    s22::AbstractVector,
    )
    n = length(X)
    length(s11) == length(s12) == length(s21) == length(s22) == n || 
        error("search_breakpoint: all vectors should have same length but " *
        "length X = $n, s11 = $(length(s11)), s12 = $(length(s12)), s21 = " *
        "$(length(s21)), s22 = $(length(s22))")

    err_optim   = typemax(Int)
    bkpts_optim = (0, 0)

    # search over all combintations of break points in two strands
    @inbounds for bkpt1 in 0:n

        # count number of errors if second haplotype is all from s22
        errors = 0
        for pos in 1:bkpt1
            if !ismissing(X[pos])
                errors += X[pos] ≠ s11[pos] + s22[pos]
            end
        end
        for pos in (bkpt1 + 1):length(X)
            if !ismissing(X[pos])
                errors += X[pos] ≠ s12[pos] +s22[pos]
            end
        end
        if errors :: Int < err_optim
            err_optim = errors
            bkpts_optim = (bkpt1, 0)

            # quick return if perfect match
            err_optim == 0 && return bkpts_optim, err_optim :: Int
        end

        # extend haplotype s21 position by position
        for bkpt2 in 1:bkpt1
            if !ismissing(X[bkpt2]) && s21[bkpt2] ≠ s22[bkpt2]
                errors -= X[bkpt2] ≠ s11[bkpt2] + s22[bkpt2]
                errors += X[bkpt2] ≠ s11[bkpt2] + s21[bkpt2]
                if errors :: Int < err_optim
                    err_optim = errors
                    bkpts_optim = (bkpt1, bkpt2)
                end
            end
        end
        for bkpt2 in (bkpt1 + 1):length(X)
            if !ismissing(X[bkpt2]) && s21[bkpt2] ≠ s22[bkpt2]
                errors -= X[bkpt2] ≠ s12[bkpt2] + s22[bkpt2]
                errors += X[bkpt2] ≠ s12[bkpt2] + s21[bkpt2]
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
