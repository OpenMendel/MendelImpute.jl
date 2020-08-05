"""
Data structure for keeping track of unique haplotypes in a window. 

- `uniqueindex`: the unique haplotype indices in a window.
- `to_unique`: map every haplotype index to the unique haplotype index
- `uniqueH`: A BitMatrix storing unique haplotypes in columns or rows, depending
    on `dims` argument of `compress_haplotypes`
- `hapmap`: Dictionary that maps every unique haplotype to all the same 
    haplotypes. If a haplotype is unique, it will not be recorded in `hapmap`
    to conserve memory.

# Example:

Suppose in window `w` there are 9 haplotypes (represented by different letters):

    [a b b c b d a d a]

then

    uniqueindex = [1, 2, 4, 6]
    to_unique = [1, 2, 2, 3, 2, 4, 1, 4, 1]
    uniqueH = [a b c d]
    hapmap = `Dict{Int32,Array{Int32,1}} with 3 entry:
        1 => [1, 7, 9]
        2 => [2, 3, 5]
        6 => [6, 8]`

where `4 => [4]` in `hapmap` has been skipped since it is a singleton.
"""
struct CompressedWindow
    uniqueindex::Vector{Int32}
    hapmap::Dict{Int32, Vector{Int32}}
    to_unique::Vector{Int32}
    uniqueH::BitMatrix
end

"""
Keeps a vector of `CompressedWindow`. Indexing off instances of 
`CompressedHaplotypes` means indexing off `CompressedHaplotypes.CW`

- `CW`: Vector of `CompressedWindow`. `CW[i]` stores unique haplotypes filtered
    with respect to all SNPs from `Hstart[i]` to `Hstart[i + 1]`
- `CW_typed`: Vector of `CompressedWindow`. `CW_typed[i]` stores unique
    haplotypes filtered with respect to typed SNPs from `Hstart[i]` to
    `Hstart[i + 1]`
- `Hstart`: `Hstart[i]` to `Hstart[i + 1]` is the range of H's SNPs (typed or
    untyped) that are in window `i`. It includes all SNP until the first 
    typed snp of window `i + 1`. 
- `X_window_range`: `X_window_range[w]` is the range of `X`'s (typed SNP's) 
    index in window `w`. 
- `sampleID`: Sample names as listed in the VCF file
- `altfreq`: Alternate allele frequency (frequency of "1" in VCF file)
- `max_unique_haplotypes`: Maximum number of unique haplotypes in each window. 
"""
struct CompressedHaplotypes
    CW::Vector{CompressedWindow}
    CW_typed::Vector{CompressedWindow}
    Hstart::Vector{Int} # we need the starting point more often than the range
    X_window_range::Vector{UnitRange}
    sampleID::Vector{String}
    chr::Vector{String}
    pos::Vector{Int}
    SNPid::Vector{Vector{String}}
    refallele::Vector{String}
    altallele::Vector{Vector{String}}
    altfreq::Vector{Float32}
    max_unique_haplotypes::Int
end
CompressedHaplotypes(windows::Int, sampleID, chr, pos, SNPid, ref, alt, altfreq,
    max_unique_haplotypes) = CompressedHaplotypes(Vector{CompressedWindow}(
    undef, windows), Vector{CompressedWindow}(undef, windows), 
    zeros(windows), [1:0 for _ in 1:windows], sampleID, chr, pos,
    SNPid, ref, alt, altfreq, max_unique_haplotypes)

nhaplotypes(x::CompressedHaplotypes) = 2length(x.sampleID)
nwindows(x::CompressedHaplotypes) = length(x.CW)
function max_width(Hunique::CompressedHaplotypes)
    win = nwindows(Hunique)
    mymax = 0
    for w in 1:win
        n = length(Hunique.X_window_range[w])
        n > mymax && (mymax = n)
    end
    return mymax
end
function max_width(reffile::String)
    endswith(reffile, ".jlso") || error("max_width can only" *
        " be called on a `CompressedHaplotypes` or `.jlso` files.")
    loaded = JLSO.load(reffile)
    compressed_Hunique = loaded[:compressed_Hunique]
    return max_width(compressed_Hunique)
end
function count_haplotypes_per_window(Hunique::CompressedHaplotypes)
    win = nwindows(Hunique)
    unique_haplotype_counts = zeros(Int, win)
    for w in 1:win
        unique_haplotype_counts[w] = length(Hunique.CW_typed[w].uniqueindex)
    end
    return unique_haplotype_counts
end
function count_haplotypes_per_window(reffile::String)
    endswith(reffile, ".jlso") || error("count_haplotypes_per_window can only" *
        " be called on a `CompressedHaplotypes` or `.jlso` files.")
    loaded = JLSO.load(reffile)
    compressed_Hunique = loaded[:compressed_Hunique]
    return count_haplotypes_per_window(compressed_Hunique)
end
avg_haplotypes_per_window(Hunique::CompressedHaplotypes) = 
    mean(count_haplotypes_per_window(Hunique))
avg_haplotypes_per_window(reffile::String) = 
    mean(count_haplotypes_per_window(reffile))

"""
For an index in unique haplotype (of typed snps), finds the first occurance of 
that haplotype in the complete reference pool for the specified window. 

This is only needed for the typed SNPs!
"""
@inline function unique_idx_to_complete_idx(unique_idx, window,
    Hunique::CompressedHaplotypes)
    return Hunique.CW_typed[window].uniqueindex[unique_idx]
end

"""
For an index in the complete haplotype pool, find its index in the unique
haplotype pool of all SNPs (typed + untyped) in specified window. 
"""
@inline function complete_idx_to_unique_all_idx(complete_idx, window, 
    Hunique::CompressedHaplotypes)
    return Hunique.CW[window].to_unique[complete_idx]
end

"""
For an index in the complete haplotype pool, find its index in the unique
haplotype pool of just the typed SNPs in specified window. 
"""
@inline function complete_idx_to_unique_typed_idx(complete_idx, window, 
    Hunique::CompressedHaplotypes)
    return Hunique.CW_typed[window].to_unique[complete_idx]
end

"""
    compress_haplotypes(reffile, tgtfile, outfile, d)

Cuts a haplotype matrix `reffile` into windows of variable width so that each
window has less than `d` unique haplotypes. Saves result to `outfile` as
a compressed binary format. All SNPs in `tgtfile` must be present in `reffile`. 

# Why is `tgtfile` required? 
The unique haplotypes in each window is computed on the typed SNPs only. 
A genotype matrix `tgtfile` is used to identify the typed SNPs. In the future, 
hopefully we can compute compressed haplotype panels for all genotyping 
platforms. 

# Inputs
* `reffile`: reference haplotype file name (ends in `.vcf` or `.vcf.gz`)
* `tgtfile`: target genotype file name (ends in `.vcf` or `.vcf.gz`)
* `outfile`: Output file name (ends in `.jlso`)
* `d`: Max number of unique haplotypes per window. 
"""
function compress_haplotypes(
    reffile::AbstractString,
    tgtfile::AbstractString,
    outfile::AbstractString,
    d::Int=1000,
    )
    endswith(outfile, ".jlso") || error("`outfile` does not end in '.jlso'")

    # import reference haplotypes
    H, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt = convert_ht(Bool, 
        reffile, trans=true, save_snp_info=true, 
        msg="importing reference data...")
    X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt = 
        VCFTools.convert_gt(UInt8, tgtfile, trans=true, save_snp_info=true, 
        msg = "Importing genotype file...")
    any(isnothing, indexin(X_pos, H_pos)) && 
        error("Found SNPs in target file that are not in reference file!")

    # compress routine
    compress_haplotypes(H, X, outfile, X_pos, H_sampleID, H_chr, H_pos, H_ids, 
        H_ref, H_alt, d)

    return nothing
end

"""
    compress_haplotypes(H, X, outfile, ...)

Compresses `H` window-by-window into `.jlso` format so that each window has `d`
unique haplotypes determined by typed SNPs position in `X`. 
"""
function compress_haplotypes(H::AbstractMatrix, X::AbstractMatrix, 
    outfile::AbstractString, X_pos::AbstractVector, H_sampleID::AbstractVector, 
    H_chr::AbstractVector, H_pos::AbstractVector, H_ids::AbstractVector, 
    H_ref::AbstractVector, H_alt::AbstractVector, d::Int)

    endswith(outfile, ".jlso") || error("`outfile` does not end in 'jlso'.")

    # some constants
    ref_snps = size(H, 1)
    tgt_snps = size(X, 1)
    Hw_idx_start = 1

    # compute window intervals based on typed SNPs
    XtoH_idx = indexin(X_pos, H_pos) # assumes all SNPs in X are in H
    Hw_typed = H[XtoH_idx, :]        # H with only typed snps
    window_ranges = get_window_intervals(Hw_typed, d)
    wins = length(window_ranges)

    # initialize compressed haplotype object
    alt_freq = reshape(sum(H, dims=2), size(H, 1)) ./ size(H, 2)
    compressed_Hunique = CompressedHaplotypes(wins, H_sampleID,
        H_chr, H_pos, H_ids, H_ref, H_alt, convert(Vector{Float32},alt_freq), d)
    compressed_Hunique.X_window_range .= window_ranges

    # record unique haplotypes and mappings window by window
    for w in 1:length(window_ranges)
        # current window ranges
        Xw_idx_start = first(window_ranges[w])
        Xw_idx_end = last(window_ranges[w])
        Xw_pos_end = X_pos[Xw_idx_end]
        Xw_pos_next_start = (w == wins ? X_pos[end] : 
            X_pos[first(window_ranges[w + 1])])
        Hw_idx_end = (w == wins ? length(H_pos) : 
            something(findnext(x -> x == Xw_pos_next_start, H_pos, 
            Hw_idx_start)) - 1)
        compressed_Hunique.Hstart[w] = Hw_idx_start

        # get current window of H
        Xw_pos = @view(X_pos[Xw_idx_start:Xw_idx_end])
        XwtoH_idx = indexin(Xw_pos, H_pos) # assumes all SNPs in X are in H
        Hw = @view(H[Hw_idx_start:Hw_idx_end, :]) # including all snps
        Hw_typed = H[XwtoH_idx, :]                # including only typed snps

        # find unique haplotypes on all SNPs
        mapping = groupslices(Hw, dims = 2)
        unique_idx = unique(mapping)
        hapmap = Dict{Int32, Vector{Int32}}()
        for idx in unique_idx
            redundant_haplotypes = findall(x -> x == idx, mapping)
            if length(redundant_haplotypes) == 1
                continue # skip singletons
            else
                hapmap[idx] = redundant_haplotypes
            end
        end
        complete_to_unique = indexin(mapping, unique_idx)
        uniqueH = Hw[:, unique_idx]
        compressed_Hunique.CW[w] = MendelImpute.CompressedWindow(unique_idx, 
            hapmap, complete_to_unique, uniqueH)

        # find unique haplotypes on typed SNPs
        mapping_typed = groupslices(Hw_typed, dims = 2)
        unique_idx_typed = unique(mapping_typed)
        hapmap_typed = Dict{Int32, Vector{Int32}}()
        for idx in unique_idx_typed
            redundant_haplotypes = findall(x -> x == idx, mapping_typed)
            if length(redundant_haplotypes) == 1
                continue # skip singletons
            else
                hapmap_typed[idx] = redundant_haplotypes
            end
        end
        complete_to_unique_typed = indexin(mapping_typed, unique_idx_typed)
        uniqueH_typed = Hw_typed[:, unique_idx_typed]
        compressed_Hunique.CW_typed[w] = 
            MendelImpute.CompressedWindow(unique_idx_typed, hapmap_typed, 
            complete_to_unique_typed, uniqueH_typed)

        # update Hw_idx_start
        Hw_idx_start = Hw_idx_end + 1
    end

    JLSO.save(outfile, :compressed_Hunique => compressed_Hunique, 
        format=:julia_serialize, compression=:gzip)

    return compressed_Hunique
end

"""
    get_window_intervals(H, d, low, high...)

Find window intervals so that every window has `d` or less unique haplotypes
by recursively dividing windows into halves. 

# Inputs
- `H`: The full haplotype matrix (on typed SNPs). Each column is a haplotype.
- `d`: Number of unique haplotypes in each window
- `low`: start of current window
- `high`: end of current window
- `intervals`: Vector of window ranges. This is also the return vector
- `seen`: storage container

# Output
- `intervals`: Vector of window ranges. This is also the return vector
"""
function get_window_intervals(
    H::AbstractMatrix,
    d::Int, 
    low::Int=1, 
    high::Int=size(H, 1),
    intervals = UnitRange[],
    seen=BitSet(),
    )
    
    unique_columns_maps = groupslices(view(H, low:high, :), dims = 2)
    k = count_unique(unique_columns_maps, seen)
    if k ≤ d
        push!(intervals, low:high)
    else
        mid = (low + high) >>> 1
        get_window_intervals(H, d, low,     mid,  intervals, seen)
        get_window_intervals(H, d, mid + 1, high, intervals, seen)
    end

    return sort!(intervals)
end

"""
    count_unique(v, seen)

Count the number of unique elements in `v`, using `seen` as storage.
"""
function count_unique(v::AbstractVector{<:Integer}, seen::AbstractSet=BitSet())
    empty!(seen)
    s = 0
    for i in v
        if i ∉ seen
            s += 1
            push!(seen, i)
        end
    end
    return s
end
