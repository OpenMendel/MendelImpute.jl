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

# """
# Data structure for keeping track of unique haplotypes in each window. 

# Let `w` be the current window, then `uniqueindex[w]` is the unique haplotype indices
# of window `w`. `hapmap1[w]` maps every haplotype to its unique haplotype. Once 
# the best haplotype pair (h1, h2) is identified, one can use `hapmap` to 
# generate a set of redundant haplotypes in window `w`. 

# # Example:

# Suppose in window `w` there are 9 haplotypes (represented by different letters):

# 	H[w, :] = [a b b c b d a d a]

# then

# 	uniqueindex[w]  = [1 2 4 6]
# 	hapmap[w] = [1 2 2 4 2 6 1 6 1].

# `hapmap` is used to find the set of matching haplotypes after identifying the best
# haplotype pair. If (b, c) is the optimal haplotype pair in the current window, then 
# the set of haplotypes that matches (b, c) is in columns {2, 3, 4, 5}
# """
# struct UniqueHaplotypeMaps
#     uniqueindex::Vector{Vector{Int}}
#     hapmap::Vector{Vector{Int}}
# end
# UniqueHaplotypeMaps(windows::Int, haps::Int) = UniqueHaplotypeMaps(Vector{Vector{Int}}(undef, windows), [zeros(Int, haps) for i in 1:windows])

"""
Data structure for holding the unique haplotypes and mappings for each window.

+ `unique_index`: Vector of BitVectors where each bitvector has length equal 
to the number of haplotypes. unique_index[w][i] = 1 if the ith haplotype in window w is unique

+ `redundant_map`: Vector of dictionaries where each dictionary stores the key/value 
mapping for redundant/unique haplotypes. redundant_map[w][i] = j means haplotype i in 
window w is not unique, and the matching unique haplotype is j. 
"""
struct UniqueHaplotypes
    unique_index::Vector{BitVector}
    redundant_map::Vector{Dict{Int64, Int64}}
end
UniqueHaplotypes(windows::Int, haps::Int) = UniqueHaplotypes([trues(haps) for i in 1:windows], [Dict{Int64, Int64}() for i in 1:windows])

# """
# Data structure for storing the redundant haplotypes matching the optimal haplotype in each window. 

# Each column is a person. Rows are `BitSet`s storing redundant haplotypes for each window 
# """
# struct RedundantHaplotypeSet
#     p::Matrix{BitSet}
# end
# RedundantHaplotypeSet(windows, people) = RedundantHaplotypeSet([BitSet() for i in 1:windows, j in 1:people])

# Base.getindex(h::RedundantHaplotypeSet, i::Int, j::Int) = h.p[i, j]
# Base.size(h::RedundantHaplotypeSet) = size(h.p)
# Base.size(h::RedundantHaplotypeSet, k::Int) = size(h.p, k)

# struct PeoplesRedundantHaplotypeSet
#     strand1::RedundantHaplotypeSet
#     strand2::RedundantHaplotypeSet
# end
# PeoplesRedundantHaplotypeSet(windows::Int, people::Int) = PeoplesRedundantHaplotypeSet(RedundantHaplotypeSet(windows, people), RedundantHaplotypeSet(windows, people))

# Base.size(h::PeoplesRedundantHaplotypeSet) = size(h.strand1)
# Base.size(h::PeoplesRedundantHaplotypeSet, k::Int) = size(h.strand1, k)

"""
Data structure for storing all haplotypes that match the optimal haplotype in each window for a person, keeping track of strand.
"""
struct OptimalHaplotypeSet
    strand1::Vector{BitVector}
    strand2::Vector{BitVector}
end

