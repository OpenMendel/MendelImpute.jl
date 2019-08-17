__precompile__()

module MendelImpute

#needed for groupslices function
import Base.hash
import Base.Cartesian, Base.Cartesian.@nloops, Base.Cartesian.@nref

using LinearAlgebra
using StatsBase
using ElasticArrays

export continue_haplotype,
    haplopair!, haplopair, haploimpute!,
    impute!, phase,
    search_breakpoint,
    unique_haplotypes, unique_haplotype_idx,
    groupslices, groupslices!, 
    phase2

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
There are a total of `length` windows. In each window, `uniqH[i]`
is an integer vector denoting the unique columns of window `i` (without
repeats). `hapmap[i]` is an integer vector where `hap_match[w][i]` 
finds the matching unique haplotype in `uniqH`. 
"""
struct UniqueHapSet
    uniqH::Vector{Vector{Int}}
    hapmap::Vector{Vector{Int}}
end
UniqueHapSet(windows, haps) = UniqueHapSet(Vector{Vector{Int}}(undef, windows), [zeros(Int, haps) for i in 1:windows])

# utilities for haplotyping
include("haplotyping.jl")

end # module
