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
    uniqueH::BitMatrix
end

"""
Keeps a vector of `CompressedWindow`. 

- `cw`: Vector of `CompressedWindow`. 
- `width`: The number of SNPs per `CompressedWindow`
"""
struct CompressedHaplotypes
    cw::Vector{CompressedWindow}
    width::Int
end
CompressedHaplotypes(windows::Int, width::Int) = CompressedHaplotypes(Vector{CompressedWindow}(undef, windows), width)

# methods to implement: https://docs.julialang.org/en/v1/manual/interfaces/#Indexing-1
Base.getindex(x::CompressedHaplotypes, w::Int) = x.cw[w]
Base.setindex!(x::CompressedHaplotypes, v::CompressedWindow, w::Int) = x.cw[w] = v
Base.firstindex(x::CompressedHaplotypes) = firstindex(x.cw)
Base.lastindex(x::CompressedHaplotypes) = lastindex(x.cw)

"""
    compress_haplotypes(vcffile, outfile, [width], [dims], [flankwidth])

For each window, finds unique haplotype indices stored in the columns of H and 
saves a mapping vector of unique columns of H. See `UniqueHaplotypeMaps` data 
structure for examples. 

# Input
* `vcffile`: file name
* `outfile`: Output file name (with or without `.jld2` extension)
* `width`: Number of SNPs per window
* `dims`: Orientation of `H`. `2` means columns of `H` are a haplotype vectors. `1` means rows of `H` are. 
* `flankwidth`: Number of SNPs flanking the sliding window (defaults to 10% of `width`)

# Output
* `hapset`: Data structure for keeping track of unique haplotypes in each window (written to `outfile`). 
"""
function compress_haplotypes(
    vcffile::AbstractString,
    outfile::AbstractString,
    width::Int;
    dims::Int, 
    flankwidth::Int = 0
    )
    # import data
    trans = (dims == 2 ? true : false)
    H = convert_ht(Bool, vcffile, trans=trans)

    # initialize constants
    dims <= 2 || error("Currently dims can only be 1 or 2, but was $dims.")
    if dims == 2
        p, d = size(H)
    elseif dims == 1
        d, p = size(H)
    end
    windows = floor(Int, p / width)

    # initialize compressed haplotype object
    hapset = CompressedHaplotypes(windows, width)

    # record unique haplotypes and mappings window by window
    for w in 1:windows
        if w == 1
            cur_range = 1:(width + flankwidth)
        elseif w == windows
            cur_range = ((windows - 1) * width - flankwidth + 1):p
        else
            cur_range = ((w - 1) * width - flankwidth + 1):(w * width + flankwidth)
        end

        H_cur_window = (dims == 2 ? view(H, cur_range, :) : view(H, :, cur_range))
        hapmap = groupslices(H_cur_window, dims=dims)
        unique_idx = unique(hapmap)
        uniqueH = convert(BitMatrix, (dims == 2 ? H_cur_window[:, unique_idx] : H_cur_window[unique_idx, :]))
        hapset[w] = CompressedWindow(unique_idx, hapmap, cur_range, uniqueH)
    end

    # save hapset to binary file using JLD2 package
    endswith(outfile, ".jld2") || (outfile = outfile * ".jld2")
    @save outfile hapset

    return hapset
end
