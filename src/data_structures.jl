# This file stores data structures used in imputation and phasing.
# For JLSO compression, the relevant data structures are in compress.jl

"""
Data structure for recording haplotype mosaic of one strand:
`start[i]` to `start[i+1]` has haplotype `haplotypelabel[i]`
in `window[i]`. The haplotype label is the column index of
`CW_typed[window[i]].uniqueH`. `start[end]` to`length` has 
haplotype `haplotypelabel[end]` in `window[end]`
"""
struct HaplotypeMosaic
    length::Int64
    start::Vector{Int64}
    window::Vector{Int32}
    haplotypelabel::Vector{Int32}
end
HaplotypeMosaic(len) = HaplotypeMosaic(len, Int64[], Int32[], Int32[])

function push_Mosaic!(x::HaplotypeMosaic, y::Tuple{Int64, T}, i) where T <: Integer
    newstart, newlabel = y[1], y[2]
    xlen = length(x.start)
    if xlen != 0
        # if i == 265
        #     println("last(x.start) = $(last(x.start)), newstart = $newstart, x.haplotypelabel[xlen] = $(x.haplotypelabel[xlen]), newlabel = $newlabel")
        # end
        # check if start occurs before previous start position, since searching
        # breakpoints for 2 consecutive windows can cause "overlaps"
        xstart = last(x.start)
        if newstart ≤ xstart
            deleteat!(x.start, xlen)
            deleteat!(x.haplotypelabel, xlen)
        end
    end
    push!(x.start, newstart)
    push!(x.haplotypelabel, newlabel)
end

function push_Mosaic!(x::HaplotypeMosaic, y::Tuple{Int64, T, T}) where T <: Integer
    newstart, newlabel, newwindow = y[1], y[2], y[3]
    xlen = length(x.start)
    if xlen != 0
        # check if start occurs before previous start position, since searching
        # breakpoints for 2 consecutive windows can cause "overlaps"
        xstart = last(x.start)
        if newstart ≤ xstart
            deleteat!(x.start, xlen)
            deleteat!(x.haplotypelabel, xlen)
            deleteat!(x.window, xlen)
        end
    end
    push!(x.start, newstart)
    push!(x.haplotypelabel, newlabel)
    push!(x.window, newwindow)
end

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
