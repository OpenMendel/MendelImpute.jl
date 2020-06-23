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
- `CWrange`: `CWrange[i]` is the range of H's SNPs that are in window `i`. It includes all SNP until the first typed snp of window `i + 1`. 
- `sampleID`: Sample names for every pair of haplotypes as listed in the VCF file
- `width`: Number of typed SNPs per window
"""
struct CompressedHaplotypes
    CW::Vector{CompressedWindow}
    CW_typed::Vector{CompressedWindow}
    CWrange::Vector{UnitRange}
    width::Int
    sampleID::Vector{String}
    chr::Vector{String}
    pos::Vector{Int}
    SNPid::Vector{Vector{String}}
    refallele::Vector{String}
    altallele::Vector{Vector{String}}
end
CompressedHaplotypes(windows::Int, width, sampleID, chr, pos, SNPid, ref, alt) = CompressedHaplotypes(Vector{CompressedWindow}(undef, windows), Vector{CompressedWindow}(undef, windows), Vector{UnitRange}(undef, windows), width, sampleID, chr, pos, SNPid, ref, alt)

nhaplotypes(x::CompressedHaplotypes) = 2length(x.sampleID)

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
    width::Int,
    )
    endswith(outfile, ".jld2") || endswith(outfile, ".jlso") || error("Unrecognized compression format: `outfile` can only end in `.jlso` or `.jld2`")

    # import reference haplotypes
    H, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt = convert_ht(Bool, reffile, trans=true, save_snp_info=true, msg="importing reference data...")
    X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt = VCFTools.convert_gt(UInt8, tgtfile, trans=true, save_snp_info=true, msg = "Importing genotype file...")
    any(isnothing, indexin(X_pos, H_pos)) && error("Found SNPs in target file that are not in reference file!")

    # compress routine
    compress_haplotypes(H, X, outfile, X_pos, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt, width)

    return nothing
end

function compress_haplotypes(H::AbstractMatrix, X::AbstractMatrix, outfile::AbstractString,
    X_pos::AbstractVector, H_sampleID::AbstractVector, H_chr::AbstractVector, 
    H_pos::AbstractVector, H_ids::AbstractVector, H_ref::AbstractVector, H_alt::AbstractVector,
    width::Int)

    endswith(outfile, ".jld2") || endswith(outfile, ".jlso") || error("Unrecognized compression format: `outfile` can only end in `.jlso` or `.jld2`")

    # some constants
    ref_snps = size(H, 1)
    tgt_snps = size(X, 1)
    windows = floor(Int, tgt_snps / width)
    Hw_idx_start = 1

    # initialize compressed haplotype object
    compressed_Hunique = MendelImpute.CompressedHaplotypes(windows, width, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt)

    # record unique haplotypes and mappings window by window
    for w in 1:windows
        # current window ranges
        Xw_idx_start = (w - 1) * width + 1
        Xw_idx_end = (w == windows ? length(X_pos) : w * width)
        Xw_pos_end = X_pos[Xw_idx_end]
        Xw_pos_next = X_pos[w * width + 1]
        Hw_idx_end = (w == windows ? length(H_pos) : 
            something(findnext(x -> x == Xw_pos_next, H_pos, Hw_idx_start)) - 1)
        compressed_Hunique.CWrange[w] = Hw_idx_start:Hw_idx_end

        # get current window of H
        Xw_pos = X_pos[Xw_idx_start:Xw_idx_end]
        XwtoH_idx = indexin(Xw_pos, H_pos) # assumes all SNPs in X are in H
        Hw = H[Hw_idx_start:Hw_idx_end, :] # including all snps
        Hw_typed = H[XwtoH_idx, :]         # including only typed snps

        # find unique haplotypes on all SNPs
        hapmap = groupslices(Hw, dims = 2)
        unique_idx = unique(hapmap)
        complete_to_unique = indexin(hapmap, unique_idx)
        uniqueH = Hw[:, unique_idx]
        compressed_Hunique.CW[w] = MendelImpute.CompressedWindow(unique_idx, hapmap, complete_to_unique, uniqueH)

        # find unique haplotypes on typed SNPs
        hapmap_typed = groupslices(Hw_typed, dims = 2)
        unique_idx_typed = unique(hapmap_typed)
        complete_to_unique_typed = indexin(hapmap_typed, unique_idx_typed)
        uniqueH_typed = Hw_typed[:, unique_idx_typed]
        compressed_Hunique.CW_typed[w] = MendelImpute.CompressedWindow(unique_idx_typed, hapmap_typed, complete_to_unique_typed, uniqueH_typed)

        # update Hw_idx_start
        Hw_idx_start = Hw_idx_end + 1
    end

    # save using JLSO or JLD2
    endswith(outfile, ".jld2") && JLD2.@save outfile compressed_Hunique
    endswith(outfile, ".jlso") && JLSO.save(outfile, :compressed_Hunique => compressed_Hunique, format=:julia_serialize, compression=:gzip)

    return nothing
end

"""
For an index in unique haplotype (of typed snps), finds the first occurance of that haplotype 
in the complete reference pool for the specified window. 

This is only needed for the typed SNPs!
"""
function unique_idx_to_complete_idx(unique_idx::Int, window::Int, Hunique::CompressedHaplotypes)
    return Hunique.CW_typed[window].uniqueindex[unique_idx]
end

"""
For an index in the complete haplotype pool, find its index in the unique haplotype pool 
of all SNPs (typed + untyped) in specified window. 
"""
function complete_idx_to_unique_all_idx(complete_idx::Int, window::Int, Hunique::CompressedHaplotypes)
    return Hunique.CW[window].to_unique[complete_idx]
end

"""
For an index in the complete haplotype pool, find its index in the unique haplotype pool 
of just the typed SNPs in specified window. 
"""
function complete_idx_to_unique_typed_idx(complete_idx::Int, window::Int, Hunique::CompressedHaplotypes)
    return Hunique.CW_typed[window].to_unique[complete_idx]
end
