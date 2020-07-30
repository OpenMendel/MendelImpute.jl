# This file stores data structures used in imputation and phasing.
# For JLSO/JLD2 compression, the relevant data structures are in compress.jl

"""
Data structure for recording haplotype mosaic of one strand:
`start[i]` to `start[i+1]` has haplotype `haplotypelabel[i]`
in `window[i]`. The haplotype label is the column index of
`CW_typed[window[i]].uniqueH`. `start[end]` to`length` has 
haplotype `haplotypelabel[end]` in `window[end]`
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
HaplotypeMosaicPair(len) = HaplotypeMosaicPair(HaplotypeMosaic(len),
    HaplotypeMosaic(len))

# reset function for dynamic programming phasing
function initialize!(x::Vector{Vector{Vector{Tuple{Int32, Int32}}}})
    n = length(x)
    win = length(x[1])
    @inbounds for i in 1:n, w in 1:win
        empty!(x[i][w])
    end
end

# resize function for dynamic programming phasing
function resize!(x::Vector{Vector{Vector{Tuple{Int32, Int32}}}}, windows::Int)
    n = length(x)
    @inbounds for i in 1:n
        Base.resize!(x[i], windows)
        sizehint!(x[i], windows)
    end
end
