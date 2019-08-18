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

Let `w` be the current window, then `uniqH[w]` is the unique haplotype indices
of window `w`. `hapmap[w]` maps every haplotype to its unique haplotype. Once 
the best haplotype pair (h1, h2) is identified, one can use `hapmap` to 
generate a set of haplotypes in window `w`. 

# Example:

Suppose in window `w` there are 9 haplotypes (represented by different letters):

	Hw = [a b b c b d a d a]

then

	uniqH[w]  = [1 2 4 6]
	hapmap[w] = [1 2 2 4 2 6 1 6 1].

If (b, c) is the optimal haplotype pair, then the set of haplotypes that matches
can (b, c) is in columns {2, 3, 4, 5}
"""
struct UniqueHaplotypeMaps
    uniqH::Vector{Vector{Int}}
    hapmap::Vector{Vector{Int}}
end
UniqueHaplotypeMaps(windows, haps) = UniqueHaplotypeMaps(Vector{Vector{Int}}(undef, windows), [zeros(Int, haps) for i in 1:windows])

