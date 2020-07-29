"""
    phase_sample!(strand1, strand2, happair1, happair2, compressed_Hunique)

For a sample, intersects redundant haplotype pairs window by window,
then search for breakpoints.

# Arguments
- `happair1`: `happair1[w]` stores best haplotype index for strand1 in window `w`
- `happair2`: `happair2[w]` stores best haplotype index for strand2 in window `w`
- `compressed_Hunique`: A `CompressedHaplotypes` object

# Optional storage argument
- `seen`: Preallocated storage container

These are needed because intersection can happen in 4 ways:
A   B      A   B
|   |  or    X
C   D      C   D
"""
function phase_sample!(
    happair1::AbstractVector{<:Integer},
    happair2::AbstractVector{<:Integer},
    compressed_Hunique::CompressedHaplotypes,
    seen::AbstractSet=BitSet(),
    )

    windows = length(strand1)
    windows == length(strand2) == length(happair1) == length(happair2) ||
         error("strand1, happair1, happair2, and happair2 have different length.")
    lifespan = (1, 1) # counter to track survival time

    # get first window's optimal haplotypes
    h1 = happair1[1]
    h2 = happair2[1]
    survivors1 = get(compressed_Hunique.CW_typed[1].hapmap, h1, h1)
    survivors2 = get(compressed_Hunique.CW_typed[1].hapmap, h2, h2)

    for w in 2:windows
        # get current window's best haplotypes
        h1 = happair1[w]
        h2 = happair2[w]
        h1set = get(compressed_Hunique.CW_typed[w].hapmap, h1, h1)
        h2set = get(compressed_Hunique.CW_typed[w].hapmap, h2, h2)

        # heuristic to decide whether cross-over is better
        crossed = false
        AC = intersect_size(survivors1, h1set)
        BD = intersect_size(survivors2, h2set)
        AD = intersect_size(survivors1, h2set)
        BC = intersect_size(survivors2, h1set)
        if AC + BD < AD + BC
            crossed = true
        end

        # update strand 1 and 2
        if crossed
            if AD == 0
                happair1[(w - lifespan[1]):(w - 1)] .= survivors1[1] # record first survivor
                lifespan[1] = 1 # reset counters
            else
                intersect!(survivors1, h1set, seen)
                lifespan[1] += 1
            end
            if BC == 0
                happair2[(w - lifespan[2]):(w - 1)] .= survivors2[1] # record first survivor
                lifespan[2] = 1 # reset counters
            else
                intersect!(survivors2, h2set, seen)
                lifespan[2] += 1
            end
        else
            if AC == 0
                happair1[(w - lifespan[1]):(w - 1)] .= survivors1[1] # record first survivor
                lifespan[1] = 1 # reset counters
            else
                intersect!(survivors1, h1set, seen)
                lifespan[1] += 1
            end
            if BD == 0
                happair2[(w - lifespan[2]):(w - 1)] .= survivors2[1] # record first survivor
                lifespan[2] = 1 # reset counters
            else
                intersect!(survivors2, h2set, seen)
                lifespan[2] += 1
            end
        end
    end

    # treat last few windows separately since intersection may not become empty
    happair1[(windows - lifespan[1]):(windows - 1)] .= survivors1[1]
    happair2[(windows - lifespan[2]):(windows - 1)] .= survivors2[1]

    return nothing
end

"""
    intersect!(v::AbstractVector, u::AbstractVector, seen::BitSet=BitSet())

Computes `v ∩ u` in place and stores result in `v`.

# Arguments
- `v`: An integer vector
- `u`: An integer vector
- `seen`: Preallocated storage container
"""
function Base.intersect!(
    v::AbstractVector{<:Integer},
    u::AbstractVector{<:Integer},
    seen::AbstractSet
    )
    empty!(seen)
    for i in u
        push!(seen, i)
    end
    for i in Iterators.reverse(eachindex(v))
        @inbounds v[i] ∉ seen && deleteat!(v, i)
    end
    nothing
end

"""
    intersect_size(v::AbstractVector, u::AbstractVector, seen::BitSet=BitSet())

Computes the size of `v ∩ u` in place. Assumes `v` is usually smaller than `u`.

# Arguments
- `v`: An integer vector
- `u`: An integer vector
- `uset`: Preallocated storage container
"""
function intersect_size(
    v::AbstractVector{<:Integer},
    u::AbstractVector{<:Integer},
    seen::AbstractSet=BitSet()
    )
    empty!(seen)
    for i in u
        push!(seen, i)
    end
    s = 0
    for i in v
        i ∈ seen && (s += 1)
    end
    return s
end

function intersect_size(
    v::AbstractVector{<:Integer},
    u::Integer,
    seen::AbstractSet=BitSet()
    )
    return u in v
end
