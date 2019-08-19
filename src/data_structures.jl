"""
Data structure for recording haplotype mosaic of one strand:
`start[i]` to `start[i+1]` has haplotype `haplotypelabel[i]`
`start[end]` to `length` has haplotype `haplotypelabel[end]`
"""
struct HaplotypeMosaic
    length::Int
    start::Vector{Int}
    haplotypelabel::Vector{Int}
end
HaplotypeMosaic(len) = HaplotypeMosaic(len, Int[], Int[])

# data structure for recording haplotype mosaic of two strands
struct HaplotypeMosaicPair
    strand1::HaplotypeMosaic
    strand2::HaplotypeMosaic
end
HaplotypeMosaicPair(len) = HaplotypeMosaicPair(HaplotypeMosaic(len), HaplotypeMosaic(len))

# needed for groupslices function
struct Prehashed
    hash::UInt
end
hash(x::Prehashed) = x.hash

"""
Data structure for keeping track of unique haplotypes in each window. 

Let `w` be the current window, then `uniqueindex[w]` is the unique haplotype indices
of window `w`. `hapmap[w]` maps every haplotype to its unique haplotype. Once 
the best haplotype pair (h1, h2) is identified, one can use `hapmap` to 
generate a set of redundant haplotypes in window `w`. 

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
end
UniqueHaplotypeMaps(windows::Int, haps::Int) = UniqueHaplotypeMaps(Vector{Vector{Int}}(undef, windows), [zeros(Int, haps) for i in 1:windows])

# Each column is a person. Rows are sets storing redundant haplotypes for each window 
struct PeoplesRedundantHaplotypeSet
    p::Matrix{Set{Int}}
end
function PeoplesRedundantHaplotypeSet(windows::Int, people::Int) 
    x = PeoplesRedundantHaplotypeSet(Matrix{Set{Int}}(undef, windows, people))
    fill!(x.p, Set{Int}())
    return x
end
Base.getindex(h::PeoplesRedundantHaplotypeSet, i::Int, j::Int) = h.p[i, j]



