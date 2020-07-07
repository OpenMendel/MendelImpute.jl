"""
Data structure for recording haplotype mosaic of one strand:
`start[i]` to `start[i+1]` has haplotype `haplotypelabel[i]`
in window `window[i]` of a `CompressedWindow`. `start[end]` to
`length` has haplotype `haplotypelabel[end]` in `window[end]`
"""
struct HaplotypeMosaic
    length::Int
    start::Vector{Int}
    window::Vector{Int}
    haplotypelabel::Vector{Int32}
end
HaplotypeMosaic(len) = HaplotypeMosaic(len, Int[], Int[], Int32[])

# data structure for recording haplotype mosaic of two strands
struct HaplotypeMosaicPair
    strand1::HaplotypeMosaic
    strand2::HaplotypeMosaic
end
HaplotypeMosaicPair(len) = HaplotypeMosaicPair(HaplotypeMosaic(len), HaplotypeMosaic(len))

"""
Data structure for keeping track of unique haplotypes in each window. 

- `uniqueindex`: the unique haplotype indices of each window.
- `hapmap`: information to map every haplotype in each window to the unique one
- `range`: range of SNPs in each window (may have overlaps in the first 2 due to flanking windows)

Let `w` be the current window, then `uniqueindex[w]` is the unique haplotype indices
of window `w`. `hapmap[w]` maps every haplotype to its unique haplotype. `range[w]`
is where the range of SNPs where the unique haplotypes are computed (may have overlaps
in beginning and end). Once the best haplotype pair (h1, h2) is identified, one can 
use `hapmap` to generate a set of redundant haplotypes in window `w`. 

# Example:

Suppose in window `w` there are 9 haplotypes (represented by different letters):

    H[w, :] = [a b b c b d a d a]

then

    uniqueindex[w]  = [1 2 4 6]
    hapmap[w] = [1 2 2 4 2 6 1 6 1].

`hapmap` is used to find the set of matching haplotypes after identifying the best
haplotype pair. If (b, c) is the optimal haplotype pair in the current window, then 
the set of haplotypes that matches (b, c) is in columns {2, 3, 4, 5}
"""
struct UniqueHaplotypeMaps
    uniqueindex::Vector{Vector{Int}}
    hapmap::Vector{Vector{Int}}
    range::Vector{UnitRange}
end
UniqueHaplotypeMaps(windows::Int, haps::Int) = UniqueHaplotypeMaps(Vector{Vector{Int}}(undef, windows), [zeros(Int, haps) for i in 1:windows], [1:0 for i in 1:windows])

"""
Data structure for storing all haplotypes that match the optimal haplotype in each window for a person, keeping track of strand.

+ `strand1[w]` stores a BitVector for window `w`. length(strand1[w]) = number of haplotypes. Entries of this BitVector is 1 if that haplotype matches the optimal haplotype.
+ `carryover1` and `carryover2` are the surviving haplotypes of the previous chunk.
"""
struct OptimalHaplotypeSet
    strand1::Vector{BitVector}
    strand2::Vector{BitVector}
    carryover1::BitVector
    carryover2::BitVector
end
OptimalHaplotypeSet(windows::Int, haps::Int) = OptimalHaplotypeSet([falses(haps) for i in 1:windows], [falses(haps) for i in 1:windows], falses(haps), falses(haps))
windows(x::OptimalHaplotypeSet) = length(x.strand1)

function initialize!(x::Vector{OptimalHaplotypeSet})
    n = length(x)
    win = windows(x[1])
    for i in 1:n
        # save last window's surviving haplotypes to carryover
        x[i].carryover1 .= x[i].strand1[win]
        x[i].carryover2 .= x[i].strand2[win]
        # reinitialize all windows to falses
        for w in 1:win
            x[i].strand1[w] .= false
            x[i].strand2[w] .= false
        end
    end
end

function resize!(x::Vector{OptimalHaplotypeSet}, windows::Int)
    n = length(x)
    for i in 1:n
        Base.resize!(x[i].strand1, windows)
        Base.resize!(x[i].strand2, windows)
        sizehint!(x[i].strand1, windows)
        sizehint!(x[i].strand2, windows)
    end
end
