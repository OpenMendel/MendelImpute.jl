"""
Data structure for keeping track of unique haplotypes in a window. 

- `uniqueindex`: the unique haplotype indices in a window.
- `hapmap`: information to map every haplotype to the first appearance of the same haplotype
- `to_unique`: map every haplotype index to the unique haplotype index
- `uniqueH`: A BitMatrix storing unique haplotypes in columns or rows, depending on `dims` argument of `compress_haplotypes`

# Example:

Suppose in window `w` there are 9 haplotypes (represented by different letters):

    [a b b c b d a d a]

then

    uniqueindex = [1 2 4 6]
    hapmap = [1 2 2 4 2 6 1 6 1]
    to_unique = [1 2 2 3 2 4 1 4 1]
    uniqueH = [a b c d]

`hapmap` is used to find the set of matching haplotypes after identifying the best
haplotype pair. If (b, c) is the optimal haplotype pair in the current window, then 
the set of haplotypes that matches (b, c) is in columns {2, 3, 4, 5}
"""
struct CompressedWindow
    uniqueindex::Vector{Int}
    hapmap::Vector{Int}
    to_unique::Vector{Int}
    uniqueH::BitMatrix
end

"""
Keeps a vector of `CompressedWindow`. Indexing off instances of `CompressedHaplotypes`
means indexing off `CompressedHaplotypes.CW`

- `CW`: Vector of `CompressedWindow`. `CW[i]` stores unique haplotypes filtered with respect to all SNPs in `CWrange[i]`
- `CW_typed`: Vector of `CompressedWindow`. `CW_typed[i]` stores unique haplotypes filtered with respect to typed SNPs in `CWrange[i]`
- `CWrange`: `CWrange[i]` is the range of SNPs that are in window `i`. 
- `sampleID`: Sample names for every pair of haplotypes as listed in the VCF file
"""
struct CompressedHaplotypes
    CW::Vector{CompressedWindow}
    CW_typed::Vector{CompressedWindow}
    CWrange::Vector{UnitRange}
    sampleID::Vector{String}
    chr::Vector{String}
    pos::Vector{Int}
    SNPid::Vector{Vector{String}}
    refallele::Vector{String}
    altallele::Vector{Vector{String}}
end
CompressedHaplotypes(windows::Int, sampleID, chr, pos, SNPid, ref, alt) = CompressedHaplotypes(Vector{CompressedWindow}(undef, windows), Vector{CompressedWindow}(undef, windows), Vector{UnitRange}(undef, windows), sampleID, chr, pos, SNPid, ref, alt)

"""
    compress_haplotypes(vcffile, tgtfile, outfile, width, [dims], [flankwidth])

For each window of `X`, finds unique haplotype indices stored in the columns of H, saves
a mapping vector to the unique col of H, and outputs compressed haplotypes as binary Julia file. 
Assumes all SNPs in `tgtfile` is present in `reffile`. 

# Inputs
* `reffile`: reference haplotype file name
* `reffile`: target genotype file name
* `outfile`: Output file name (ends in `.jld2` or `.jlso` (recommended))
* `width`: Number of typed SNPs per window. Number of SNPs in last window may be in `[width, 2width]`.

# Optional inputs
* `flankwidth`: Number of SNPs flanking the sliding window (defaults to 10% of `width`)
"""
function compress_haplotypes(
    reffile::AbstractString,
    tgtfile::AbstractString,
    outfile::AbstractString,
    width::Int
    )
    endswith(outfile, ".jld2") || endswith(outfile, ".jlso") || error("Unrecognized compression format: `outfile` can only end in `.jlso` or `.jld2`")

    # import data
    X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt = VCFTools.convert_gt(UInt8, tgtfile, trans=true, save_snp_info=true, msg = "Importing genotype file...")
    H, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt = convert_ht(Bool, reffile, trans=true, save_snp_info=true, msg="importing vcf data...")
    any(isnothing, indexin(X_pos, H_pos)) && error("Found SNPs in target file that are not in reference file!")

    # some constants
    ref_snps = size(H, 1)
    tgt_snps = size(X, 1)
    windows = floor(Int, tgt_snps / width)
    Hw_idx_start = 1

    # initialize compressed haplotype object
    compressed_Hunique = MendelImpute.CompressedHaplotypes(windows, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt)

    # record unique haplotypes and mappings window by window
    for w in 1:windows
        # current window ranges
        Xw_idx_start = (w - 1) * width + 1
        Xw_idx_end = (w == windows ? length(X_pos) : w * width)
        Xw_pos_end = X_pos[Xw_idx_end]
        Hw_idx_end = (w == windows ? length(H_pos) : something(findnext(x -> x == Xw_pos_end, H_pos, Hw_idx_start)))

        # get current window of H, including all snps and only typed snps
        Xw_pos = X_pos[Xw_idx_start:Xw_idx_end]
        XwtoH_idx = indexin(Xw_pos, H_pos)
        Hw = H[Hw_idx_start:Hw_idx_end, :] # all snps
        Hw_typed = H[XwtoH_idx, :]         # only typed snps

        # find unique haplotypes on all SNPs
        hapmap = groupslices(Hw, dims = 2)
        unique_idx = unique(hapmap)
        complete_to_unique = indexin(hapmap, unique_idx)
        uniqueH = Hw[:, unique_idx]
        compressed_Hunique.CW[w] = MendelImpute.CompressedWindow(unique_idx, hapmap, complete_to_unique, uniqueH)

        # find unique haplotypes on typed SNPs
        hapmap = groupslices(Hw_typed, dims = 2)
        unique_idx = unique(hapmap)
        complete_to_unique = indexin(hapmap, unique_idx)
        uniqueH = Hw_typed[:, unique_idx]
        compressed_Hunique.CW_typed[w] = MendelImpute.CompressedWindow(unique_idx, hapmap, complete_to_unique, uniqueH)

        # update Hw_idx_start
        Hw_idx_start = Hw_idx_end + 1
    end

    # save using JLSO or JLD2
    endswith(outfile, ".jld2") && JLD2.@save outfile compressed_Hunique
    endswith(outfile, ".jlso") && JLSO.save(outfile, :compressed_Hunique => compressed_Hunique, format=:julia_serialize, compression=:gzip)

    return compressed_Hunique
end

"""
For an index in unique haplotype, finds the first occurance of that haplotype 
in the complete reference pool for the specified window.
"""
function unique_idx_to_complete_idx(unique_idx::Int, window::Int, Hunique::CompressedHaplotypes)
    return Hunique[window].uniqueindex[unique_idx]
end

"""
For an index in the complete haplotype pool, find its index in the unique haplotype pool 
in specified window. 
"""
function complete_idx_to_unique_idx(complete_idx::Int, window::Int, Hunique::CompressedHaplotypes)
    return Hunique[window].to_unique[complete_idx]
end
