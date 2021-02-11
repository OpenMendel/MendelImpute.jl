###### This file is part of the MendelImpute.jl package.
###### It contains relevant code to compress a reference haplotype panel into
###### a `.jlso` or `.jld2` compressed panel. 

"""
Data structure for keeping track of unique haplotypes in a window. 

# Every window contains:
- `uniqueindex`: the unique haplotype indices
- `to_unique`: vertor that maps every haplotype index to the corresponding 
    unique haplotype index
- `uniqueH`: A BitMatrix storing unique haplotypes in columns
- `hapmap`: Dictionary that maps every unique haplotype to its equivalent
    haplotypes. Singleton haplotypes will not be recorded to conserve memory.

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
A compressed version of the original reference haplotype panel. Internally
it is just 2 vectors of `CompressedWindow`, some pointers for coordination,
and basic sample information.

- `CW`: Vector of `CompressedWindow`. `CW[i]` stores unique haplotypes filtered
    with respect to all SNPs from `Hstart[i]` to `Hstart[i + 1]`
- `CW_typed`: Vector of `CompressedWindow`. `CW_typed[i]` stores unique
    haplotypes filtered with respect to typed SNPs from `Hstart[i]` to
    `Hstart[i + 1]`
- `Hstart`: `Hstart[i]` to `Hstart[i + 1]` is the range of H's SNPs (typed or
    untyped) that are in window `i`. It includes all SNP until the first 
    typed snp of window `i + 1`. Used to directly index into `H` like
    `H[Hstart[i]:Hstart[i+1], :]`
- `X_window_range`: `idx = X_window_range[w]` is the range of `X`'s (typed
    SNP's) index in window `w`. Used to directly index into `X` like `X[idx, :]`.
    Does not include overlapping SNPs. 
- `sampleID`: Sample names as listed in the VCF file
- `altfreq`: Alternate allele frequency (frequency of "1" in VCF file)
- `max_unique_haplotypes`: Maximum number of unique haplotypes in each window. 
- `overlap`: Boolean indicating whether adjacent windows of typed SNPs overlap. 
"""
struct CompressedHaplotypes
    CW::Vector{CompressedWindow}
    CW_typed::Vector{CompressedWindow}
    Hstart::Vector{Int} # we need the starting point more often than the range
    X_window_range::Vector{UnitRange{Int}}
    sampleID::Vector{String}
    chr::Vector{String}
    pos::Vector{Int}
    SNPid::Vector{Vector{String}}
    refallele::Vector{String}
    altallele::Vector{Vector{String}}
    altfreq::Vector{Float32}
    max_unique_haplotypes::Int
    overlap::Bool
end
CompressedHaplotypes(windows::Int, sampleID, chr, pos, SNPid, ref, alt, altfreq,
    max_unique_haplotypes, overlap) = CompressedHaplotypes(
    Vector{CompressedWindow}( undef, windows), Vector{CompressedWindow}(undef,
    windows), zeros(windows), [1:0 for _ in 1:windows], sampleID, chr, pos,
    SNPid, ref, alt, altfreq, max_unique_haplotypes, overlap)

nhaplotypes(x::CompressedHaplotypes) = 2length(x.sampleID)
nwindows(x::CompressedHaplotypes) = length(x.CW)

"""
    max_dim(reffile::String)

Searches each compressed window of the typed SNPs and outputs the maximum window
width (i.e. number of SNPs) and maximum unique haplotypes per window. If windows
are overlapping, extra padding will be automatically inserted
"""
function max_dim(reffile::String)
    endswith(reffile, ".jlso") || error("max_width can only" *
        " be called on a `CompressedHaplotypes` or `.jlso` files.")
    loaded = JLSO.load(reffile)
    compressed_Hunique = loaded[:compressed_Hunique]
    return max_dim(compressed_Hunique)
end
function max_dim(Hunique::CompressedHaplotypes)
    win = nwindows(Hunique)
    maxwidth = 0
    max_d = 0
    for w in 1:win
        n = size(Hunique.CW_typed[w].uniqueH, 1)
        n > maxwidth && (maxwidth = n)

        d = size(Hunique.CW_typed[w].uniqueH, 2)
        d > max_d && (max_d = d)
    end
    return maxwidth, max_d
end

"""
    get_window_widths(reffile::String)

Searches each compressed window of the typed SNPs and output their window widths
(i.e. number of SNPs per window) 
"""
function get_window_widths(reffile::String)
    endswith(reffile, ".jlso") || error("max_width can only" *
        " be called on a `CompressedHaplotypes` or `.jlso` files.")
    loaded = JLSO.load(reffile)
    compressed_Hunique = loaded[:compressed_Hunique]
    return get_window_widths(compressed_Hunique)
end
function get_window_widths(Hunique::CompressedHaplotypes)
    win = nwindows(Hunique)
    widths = zeros(Int, win)
    for w in 1:win
        widths[w] = length(Hunique.X_window_range[w])
    end
    return widths
end

"""
    count_haplotypes_per_window(Hunique::CompressedHaplotypes)

Searches each compressed window of the typed SNPs and output number of unique
haplotypes per window
"""
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

"""
    avg_haplotypes_per_window(Hunique::CompressedHaplotypes)

Searches each compressed window of the typed SNPs and output average number of
unique haplotypes per window
"""
avg_haplotypes_per_window(Hunique::CompressedHaplotypes) = 
    mean(count_haplotypes_per_window(Hunique))
avg_haplotypes_per_window(reffile::String) = 
    mean(count_haplotypes_per_window(reffile))

"""
For an index in unique haplotype (of typed snps), finds the first occurance of 
that haplotype in the complete reference pool for the specified window. 
"""
function unique_idx_to_complete_idx(unique_idx, window,
    Hunique::CompressedHaplotypes)
    return Hunique.CW_typed[window].uniqueindex[unique_idx]
end

"""
For an index in unique haplotype (of all snps), finds the first occurance of 
that haplotype in the complete reference pool for the specified window. 
"""
function unique_all_idx_to_complete_idx(unique_idx, window,
    Hunique::CompressedHaplotypes)
    return Hunique.CW[window].uniqueindex[unique_idx]
end

"""
For an index in unique haplotype (of all snps), finds the haplotype in the 
unique haplotype (of typed SNPs only) pool for the specified window. 
"""
function unique_all_idx_to_unique_typed_idx(all_idx, window,
    Hunique::CompressedHaplotypes)
    complete_idx = unique_all_idx_to_complete_idx(all_idx, window, Hunique)
    return complete_idx_to_unique_typed_idx(complete_idx, window, Hunique)
end

"""
For an index in the complete haplotype pool, find its index in the unique
haplotype pool of all SNPs (typed + untyped) in specified window. 
"""
function complete_idx_to_unique_all_idx(complete_idx, window, 
    Hunique::CompressedHaplotypes)
    return Hunique.CW[window].to_unique[complete_idx]
end

"""
For an index in the complete haplotype pool, find its index in the unique
haplotype pool of just the typed SNPs in specified window. 
"""
function complete_idx_to_unique_typed_idx(complete_idx, window, 
    Hunique::CompressedHaplotypes)
    return Hunique.CW_typed[window].to_unique[complete_idx]
end

"""
    extend_to_overlap_range(Hunique::CompressedHaplotypes, window::Int, overlap::Bool)

If genotype window overlap, the unique haplotype matrix will be larger than the
genotype matrix in window `w`. This function computes the range of genotype
matrix so that the ranges match up.  
"""
function extend_to_overlap_range(
    Hunique::CompressedHaplotypes, 
    window::Int, 
    overlap::Bool
    )
    winranges = Hunique.X_window_range
    Hw = Hunique.CW_typed[window].uniqueH
    Xrange = winranges[window]
    if overlap
        extra_width = size(Hw, 1) - length(winranges[window])
        if window == 1 # first window
            Xstart = first(winranges[window])
            Xend = last(winranges[window]) + extra_width
        elseif window == length(winranges) # last window
            Xstart = first(winranges[window]) - extra_width
            Xend = last(winranges[window])
        else
            Xstart = first(winranges[window]) - (extra_width >> 1)
            Xend = last(winranges[window]) + (extra_width >> 1)
        end
        Xrange = Xstart:Xend
    end
    return Xrange
end

"""
    nonoverlap_range(Hunique::CompressedHaplotypes, window::Int, overlap::Bool)

If windows overlap, the unique haplotype matrix will include extra SNPs on 
both ends as paddings. This function calculates the range of SNPs of the unique
haplotype matrix that are not in the overlapping regions. 
"""
function nonoverlap_range(
    Hunique::CompressedHaplotypes, 
    window::Int
    )
    winranges = Hunique.X_window_range
    Hw = Hunique.CW_typed[window].uniqueH
    Hw_width = size(Hw, 1)
    extra_width = Hw_width - length(winranges[window])

    if window == 1 # first window
        Hstart = 1
        Hend = Hw_width - extra_width
    elseif window == length(winranges) # last window
        Hstart = 1 + extra_width
        Hend = Hw_width
    else
        Hstart = 1 + (extra_width >> 1)
        Hend = Hw_width - (extra_width >> 1)
    end
    
    return Hstart:Hend
end

"""
    compress_haplotypes(reffile::String, tgtfile::String, outfile::String, 
        [d::Int], [minwidth::Int], [overlap::Float64])

Cuts a haplotype matrix `reffile` into windows of variable width so that each
window has less than `d` unique haplotypes. Saves result to `outfile` as
a compressed binary format. All SNPs in `tgtfile` must be present in `reffile`. 

# Why is `tgtfile` required? 
The unique haplotypes in each window is computed on the typed SNPs only. 
A genotype matrix `tgtfile` is used to identify the typed SNPs. In the future, 
hopefully we can pre-compute compressed haplotype panels for all genotyping 
platforms and provide them as downloadable files. But currently, users must
run this function by themselves. 

# Inputs
* `reffile`: reference haplotype file name (ends in `.vcf` or `.vcf.gz`)
* `tgtfile`: VCF or PLINK files. VCF files should end in `.vcf` or `.vcf.gz`.
    PLINK files should exclude `.bim/.bed/.fam` suffixes but the trio must all
    be present in the same directory.
* `outfile`: Output file name (ends in `.jlso`)

# Optional Inputs
* `d`: Max number of unique haplotypes per genotype window (default `d = 1000`). 
* `minwidth`: Minimum number of typed SNPs per window (default 0)
* `overlap`: How much overlap between adjacent genotype windows in percentage of
    each window's width (default 0.0)
"""
function compress_haplotypes(
    reffile::AbstractString,
    tgtfile::AbstractString,
    outfile::AbstractString,
    d::Int=1000,
    minwidth::Int=0,
    overlap::Float64=0.0
    )
    # first handle error
    endswith(outfile, ".jld2") || endswith(outfile, ".jlso") || 
        error("Unrecognized compression format: `outfile` can only end in " * 
        "`.jlso` or `.jld2`")
    0.0 ≤ overlap ≤ 1.0 || error("overlap must be a percentage")

    if endswith(tgtfile, ".vcf") || endswith(tgtfile, ".vcf.gz")
        X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt = 
            VCFTools.convert_gt(UInt8, tgtfile, trans=true, 
            save_snp_info=true, msg = "Importing genotype file...")
    elseif isplink(tgtfile)
        # convert SnpArray data to matrix.
        X_snpdata = SnpArrays.SnpData(tgtfile)
        X = convert(Matrix{Union{UInt8, Missing}}, X_snpdata.snparray')
        X[findall(isone, X)] .= missing     # 0x01 encodes missing
        X[findall(x -> x === 0x02, X)] .= 1 # 0x02 is 1
        X[findall(x -> x === 0x03, X)] .= 2 # 0x03 is 2
        # get other relevant information
        X_sampleID = X_snpdata.person_info[!, :iid]
        X_chr = X_snpdata.snp_info[!, :chromosome]
        X_pos = X_snpdata.snp_info[!, :position]
        X_ids = X_snpdata.snp_info[!, :snpid]
        X_ref = X_snpdata.snp_info[!, :allele1]
        X_alt = X_snpdata.snp_info[!, :allele2]
    else
        error("Unrecognized target file format: target file can only be VCF" *
            " files (ends in .vcf or .vcf.gz) or PLINK files (do not include" *
            " .bim/bed/fam and all three files must exist in 1 directory)")
    end    
    
    # import reference haplotypes
    H, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt = convert_ht(Bool, 
        reffile, trans=true, save_snp_info=true, 
        msg="importing reference data...")
    any(isnothing, indexin(X_pos, H_pos)) && error("Found SNPs in target " * 
        "file that are not in reference file! Please filter them out first!")

    # compress routine
    compress_haplotypes(H, outfile, X_pos, H_sampleID, H_chr, H_pos, H_ids, 
        H_ref, H_alt, d, minwidth, overlap)

    return nothing
end

function compress_haplotypes(H::AbstractMatrix, outfile::AbstractString, 
    X_pos::AbstractVector, H_sampleID::AbstractVector, H_chr::AbstractVector,
    H_pos::AbstractVector, H_ids::AbstractVector, H_ref::AbstractVector,
    H_alt::AbstractVector, d::Int, minwidth::Int, overlap::Float64)

    endswith(outfile, ".jld2") || endswith(outfile, ".jlso") || 
        error("Unrecognized compression format: `outfile` can only end in " * 
        "`.jlso` or `.jld2`")
    0.0 ≤ overlap ≤ 1.0 || error("overlap must be a percentage")

    # some constants
    ref_snps = size(H, 1)
    tgt_snps = length(X_pos)

    # compute window intervals based on typed SNPs
    XtoH_idx = indexin(X_pos, H_pos) # assumes all SNPs in X are in H
    any(isnothing, XtoH_idx) && error("Found SNPs in target " * 
        "file that are not in reference file! Please filter them out first!")
    Hw_typed = H[XtoH_idx, :]        # H with only typed snps
    window_ranges, mapping_typed = get_window_intervals(Hw_typed, d, minwidth)
    wins = length(window_ranges)

    # initialize compressed haplotype object
    alt_freq = reshape(sum(H, dims=2), size(H, 1)) ./ size(H, 2)
    compressed_Hunique = CompressedHaplotypes(wins, H_sampleID,
        H_chr, H_pos, H_ids, H_ref, H_alt, convert(Vector{Float32},alt_freq),
        d, overlap > 0.0)
    compressed_Hunique.X_window_range .= window_ranges

    # record unique haplotypes and mappings window by window
    Hw_idx_start = 1
    for w in 1:length(window_ranges)
        # genotype matrix's current window ranges
        Xw_idx_start = first(window_ranges[w])
        Xw_idx_end = last(window_ranges[w])
        if overlap > 0.0 # overlaps between adjacent windows
            flankwidth = round(Int, overlap*(Xw_idx_end-Xw_idx_start), RoundUp)
            w == 1 || (Xw_idx_start -= flankwidth)
            w == wins || (Xw_idx_end += flankwidth)
        end

        # haplotype reference panel's window ranges
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
        Hw_typed = @view(H[XwtoH_idx, :])         # including only typed snps

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
        unique_idx_typed = unique(mapping_typed[w])
        hapmap_typed = Dict{Int32, Vector{Int32}}()
        for idx in unique_idx_typed
            redundant_haplotypes = findall(x -> x == idx, mapping_typed[w])
            if length(redundant_haplotypes) == 1
                continue # skip singletons
            else
                hapmap_typed[idx] = redundant_haplotypes
            end
        end
        complete_to_unique_typed = indexin(mapping_typed[w], unique_idx_typed)
        uniqueH_typed = Hw_typed[:, unique_idx_typed]
        compressed_Hunique.CW_typed[w] = 
            MendelImpute.CompressedWindow(unique_idx_typed, hapmap_typed, 
            complete_to_unique_typed, uniqueH_typed)

        # update Hw_idx_start
        Hw_idx_start = Hw_idx_end + 1
    end

    # save using JLSO or JLD2
    if endswith(outfile, ".jld2")
        jldopen(outfile, "w", compress=true) do file
            file["compressed_Hunique"] = compressed_Hunique
        end
    else #jlso
        JLSO.save(outfile, :compressed_Hunique => 
            compressed_Hunique, format=:julia_serialize, compression=:gzip)
    end

    return compressed_Hunique
end

"""
    get_window_intervals(H, d, minwidth, low, high, intervals, seen)

Find window intervals so that every window has `d` or less unique haplotypes
by recursively dividing windows into halves. 

# Inputs
- `H`: The full haplotype matrix (on typed SNPs). Each column is a haplotype.
- `d`: Number of unique haplotypes in each window
- `minwidth`: Minimum number of SNPs per window (default 20)
- `low`: start of current window
- `high`: end of current window
- `intervals`: Vector of window ranges. This is also the return vector
- `seen`: storage container

# Output
- `intervals`: Vector of window ranges. Each range have ≤ `d` unique haplotypes
"""
function get_window_intervals(
    H::AbstractMatrix,
    d::Int, 
    minwidth::Int=20,
    low::Int=1, 
    high::Int=size(H, 1),
    intervals = UnitRange[],
    unique_columns_maps = Vector{Vector{Int}}(),
    seen=BitSet(),
    )
    mid = (low + high) >> 1
    unique_columns_map = groupslices(view(H, low:high, :), dims = 2)
    k = count_unique(unique_columns_map, seen)
    if k ≤ d || length(low:mid) ≤ minwidth
        push!(intervals, low:high)
        push!(unique_columns_maps, unique_columns_map)
    else
        get_window_intervals(H, d, minwidth, low,     mid,  intervals, unique_columns_maps, seen)
        get_window_intervals(H, d, minwidth, mid + 1, high, intervals, unique_columns_maps, seen)
    end

    return intervals, unique_columns_maps
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

"""
    read_jlso(file::AbstractString)

Imports a `.jlso`-compressed reference haplotype panel.
"""
function read_jlso(reffile::AbstractString)
    loaded = JLSO.load(reffile)
    return loaded[:compressed_Hunique]
end
