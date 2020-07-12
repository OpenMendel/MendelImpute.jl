"""
Data structure for keeping track of unique haplotypes in a window. 

- `uniqueindex`: the unique haplotype indices in a window.
- `to_unique`: map every haplotype index to the unique haplotype index
- `uniqueH`: A BitMatrix storing unique haplotypes in columns or rows, depending on `dims` argument of `compress_haplotypes`
- `hapmap`: Dictionary that maps every unique haplotype to all the same haplotypes. If a haplotype is unique, it will not be recorded in `hapmap` to conserve memory.

# Example:

Suppose in window `w` there are 9 haplotypes (represented by different letters):

    [a b b c b d a d a]

then

    uniqueindex = [1, 2, 4, 6]
    to_unique = [1, 2, 2, 3, 2, 4, 1, 4, 1]
    uniqueH = [a b c d]
    hapmap = `Dict{Int64,Array{Int64,1}} with 3 entry:
        1 => [1, 7, 9]
        2 => [2, 3, 5]
        6 => [6, 8]`

where `4 => [4]` in `hapmap` has been skipped since it is a singleton.
"""
struct CompressedWindow
    uniqueindex::Vector{Int}
    hapmap::Dict{Int64, Vector{Int64}}
    to_unique::Vector{Int}
    uniqueH::BitMatrix
end

"""
Keeps a vector of `CompressedWindow`. Indexing off instances of `CompressedHaplotypes`
means indexing off `CompressedHaplotypes.CW`

- `CW`: Vector of `CompressedWindow`. `CW[i]` stores unique haplotypes filtered with respect to all SNPs from `start[i]` to `start[i + 1]`
- `CW_typed`: Vector of `CompressedWindow`. `CW_typed[i]` stores unique haplotypes filtered with respect to typed SNPs from `start[i]` to `start[i + 1]`
- `start`: `start[i]` to `start[i + 1]` is the range of H's SNPs that are in window `i`. It includes all SNP until the first typed snp of window `i + 1`. 
- `sampleID`: Sample names for every pair of haplotypes as listed in the VCF file
- `width`: Number of typed SNPs per window
- `altfreq`: Alternate allele frequency (frequency of "1" in VCF file)
"""
struct CompressedHaplotypes
    CW::Vector{CompressedWindow}
    CW_typed::Vector{CompressedWindow}
    start::Vector{Int}
    width::Int
    sampleID::Vector{String}
    chr::Vector{String}
    pos::Vector{Int}
    SNPid::Vector{Vector{String}}
    refallele::Vector{String}
    altallele::Vector{Vector{String}}
    altfreq::Vector{Float32}
end
CompressedHaplotypes(windows::Int, width, sampleID, chr, pos, SNPid, ref, alt, altfreq) = CompressedHaplotypes(Vector{CompressedWindow}(undef, windows), Vector{CompressedWindow}(undef, windows), zeros(windows), width, sampleID, chr, pos, SNPid, ref, alt, altfreq)

nhaplotypes(x::CompressedHaplotypes) = 2length(x.sampleID)
windows(x::CompressedHaplotypes) = length(x.CW)
function count_haplotypes_per_window(Hunique::CompressedHaplotypes)
    win = windows(Hunique)
    unique_haplotype_counts = zeros(Int, win)
    for w in 1:win
        unique_haplotype_counts[w] = length(Hunique.CW_typed[w].uniqueindex)
    end
    return unique_haplotype_counts
end
function count_haplotypes_per_window(reffile::String)
    endswith(reffile, ".jlso") || error("count_haplotypes_per_window can only be called on a `CompressedHaplotypes` or `.jlso` files.")
    loaded = JLSO.load(reffile)
    compressed_Hunique = loaded[:compressed_Hunique]
    return count_haplotypes_per_window(compressed_Hunique)
end
avg_haplotypes_per_window(Hunique::CompressedHaplotypes) = mean(count_haplotypes_per_window(Hunique))
avg_haplotypes_per_window(reffile::String) = mean(count_haplotypes_per_window(reffile))


"""
For an index in unique haplotype (of typed snps), finds the first occurance of that haplotype 
in the complete reference pool for the specified window. 

This is only needed for the typed SNPs!
"""
@inline function unique_idx_to_complete_idx(unique_idx, window, Hunique::CompressedHaplotypes)
    return Hunique.CW_typed[window].uniqueindex[unique_idx]
end

"""
For an index in the complete haplotype pool, find its index in the unique haplotype pool 
of all SNPs (typed + untyped) in specified window. 
"""
@inline function complete_idx_to_unique_all_idx(complete_idx, window, Hunique::CompressedHaplotypes)
    return Hunique.CW[window].to_unique[complete_idx]
end

"""
For an index in the complete haplotype pool, find its index in the unique haplotype pool 
of just the typed SNPs in specified window. 
"""
@inline function complete_idx_to_unique_typed_idx(complete_idx, window, Hunique::CompressedHaplotypes)
    return Hunique.CW_typed[window].to_unique[complete_idx]
end

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

"""
`X` and `H` stores genotypes/haplotypes in columns. 
"""
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
    alt_freq = reshape(sum(H, dims=2), size(H, 1)) ./ size(H, 2)
    compressed_Hunique = MendelImpute.CompressedHaplotypes(windows, width, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt, convert(Vector{Float32}, alt_freq))

    # record unique haplotypes and mappings window by window
    for w in 1:windows
        # current window ranges
        Xw_idx_start = (w - 1) * width + 1
        Xw_idx_end = (w == windows ? length(X_pos) : w * width)
        Xw_pos_end = X_pos[Xw_idx_end]
        Xw_pos_next_start = (w == windows ? X_pos[end] : X_pos[w * width + 1])
        Hw_idx_end = (w == windows ? length(H_pos) : 
            something(findnext(x -> x == Xw_pos_next_start, H_pos, Hw_idx_start)) - 1)
        compressed_Hunique.start[w] = Hw_idx_start

        # get current window of H
        Xw_pos = X_pos[Xw_idx_start:Xw_idx_end]
        XwtoH_idx = indexin(Xw_pos, H_pos) # assumes all SNPs in X are in H
        Hw = H[Hw_idx_start:Hw_idx_end, :] # including all snps
        Hw_typed = H[XwtoH_idx, :]         # including only typed snps

        # find unique haplotypes on all SNPs
        mapping = groupslices(Hw, dims = 2)
        unique_idx = unique(mapping)
        hapmap = Dict{Int, Vector{Int}}()
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
        compressed_Hunique.CW[w] = MendelImpute.CompressedWindow(unique_idx, hapmap, complete_to_unique, uniqueH)

        # find unique haplotypes on typed SNPs
        mapping_typed = groupslices(Hw_typed, dims = 2)
        unique_idx_typed = unique(mapping_typed)
        hapmap_typed = Dict{Int, Vector{Int}}()
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
        compressed_Hunique.CW_typed[w] = MendelImpute.CompressedWindow(unique_idx_typed, hapmap_typed, complete_to_unique_typed, uniqueH_typed)

        # update Hw_idx_start
        Hw_idx_start = Hw_idx_end + 1
    end

    # save using JLSO or JLD2
    endswith(outfile, ".jld2") && JLD2.@save outfile compressed_Hunique
    endswith(outfile, ".jlso") && JLSO.save(outfile, :compressed_Hunique => compressed_Hunique, format=:julia_serialize, compression=:gzip)

    return nothing
end
