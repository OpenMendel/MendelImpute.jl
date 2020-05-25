"""
    Data structure saving each VCF record's relevant information
"""
struct VCFInfo
    chr::Vector{String}
    pos::Vector{Int}
    SNPid::Vector{Vector{String}}
    refallele::Vector{String}
    altallele::Vector{Vector{String}}
end

"""
Data structure for keeping track of unique haplotypes in a window. 

- `uniqueindex`: the unique haplotype indices in a window.
- `hapmap`: information to map every haplotype to the first appearance of the same haplotype
- `range`: range of SNPs of the window
- `uniqueH`: A BitMatrix storing unique haplotypes in columns or rows, depending on `dims` argument of `compress_haplotypes`

# Example:

Suppose in window `w` there are 9 haplotypes (represented by different letters):

    [a b b c b d a d a]

then

    uniqueindex = [1 2 4 6]
    hapmap = [1 2 2 4 2 6 1 6 1]
    uniqueH = [a b c d]

`hapmap` is used to find the set of matching haplotypes after identifying the best
haplotype pair. If (b, c) is the optimal haplotype pair in the current window, then 
the set of haplotypes that matches (b, c) is in columns {2, 3, 4, 5}
"""
struct CompressedWindow
    uniqueindex::Vector{Int}
    hapmap::Vector{Int}
    range::UnitRange
    vcfinfo::VCFInfo
    uniqueH::BitMatrix
end

"""
Keeps a vector of `CompressedWindow`. 

- `cw`: Vector of `CompressedWindow`. 
- `width`: The number of SNPs per `CompressedWindow`. The last window may have a different width
- `snps`: Total number of snps in every window
- `sampleID`: Sample names for every pair of haplotypes as listed in the VCF file
"""
struct CompressedHaplotypes
    cw::Vector{CompressedWindow}
    width::Int
    snps::Int
    sampleID::Vector{String}
end
CompressedHaplotypes(windows::Int, width::Int, snps::Int, sampleID::Vector{String}) = CompressedHaplotypes(Vector{CompressedWindow}(undef, windows), width, snps, sampleID)

# methods to implement: https://docs.julialang.org/en/v1/manual/interfaces/#Indexing-1
Base.getindex(x::CompressedHaplotypes, w::Int) = x.cw[w]
Base.setindex!(x::CompressedHaplotypes, v::CompressedWindow, w::Int) = x.cw[w] = v
Base.firstindex(x::CompressedHaplotypes) = firstindex(x.cw)
Base.lastindex(x::CompressedHaplotypes) = lastindex(x.cw)

"""
    compress_haplotypes(vcffile, outfile, [width], [dims], [flankwidth])

For each window, finds unique haplotype indices stored in the columns/rows of H, saves
a mapping vector to the unique col/row of H, and outputs result as binary `.jld2` file. 

# Input
* `vcffile`: file name
* `outfile`: Output file name (with or without `.jld2` extension)
* `width`: Number of SNPs per window
* `dims`: Orientation of `H`. `2` means columns of `H` are a haplotype vectors. `1` means rows of `H` are. 
* `flankwidth`: Number of SNPs flanking the sliding window (defaults to 10% of `width`)
"""
function compress_haplotypes(
    vcffile::AbstractString,
    outfile::AbstractString,
    width::Int;
    dims::Int, 
    flankwidth::Int = 0
    )
    # import data
    trans = (dims == 2 ? true : error("currently VCFTools only import chr/pos...info when H transposed"))
    H, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt = convert_ht(Bool, vcffile, trans=trans, save_snp_info=true, msg="importing vcf data...")

    # some constants
    snps = (dims == 2 ? size(H, 1) : size(H, 2))
    windows = floor(Int, snps / width)

    # initialize compressed haplotype object
    compressed_Hunique = CompressedHaplotypes(windows, width, snps, H_sampleID)

    # record unique haplotypes and mappings window by window
    for w in 1:windows
        if w == 1
            cur_range = 1:(width + flankwidth)
        elseif w == windows
            cur_range = ((windows - 1) * width - flankwidth + 1):snps
        else
            cur_range = ((w - 1) * width - flankwidth + 1):(w * width + flankwidth)
        end

        H_cur_window = (dims == 2 ? view(H, cur_range, :) : view(H, :, cur_range))
        hapmap = groupslices(H_cur_window, dims=dims)
        unique_idx = unique(hapmap)
        uniqueH = (dims == 2 ? H_cur_window[:, unique_idx] : H_cur_window[unique_idx, :])
        info = VCFInfo(H_chr[cur_range], H_pos[cur_range], H_ids[cur_range], H_ref[cur_range], H_alt[cur_range])
        compressed_Hunique[w] = CompressedWindow(unique_idx, hapmap, cur_range, info, uniqueH)
    end

    # save to binary file using JLD2 package
    endswith(outfile, ".jld2") || (outfile = outfile * ".jld2")
    @save outfile compressed_Hunique

    return compressed_Hunique
end
