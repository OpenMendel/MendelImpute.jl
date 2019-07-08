__precompile__()

module MendelImpute

using LinearAlgebra

export continue_haplotype,
    haplopair!, haplopair, haploimpute!,
    impute!, phase,
    search_breakpoint,
    filter_redundant_haplotypes, concats, redundant_index

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

# utilities for haplotyping
include("haplotyping.jl")

end # module
