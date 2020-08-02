"""
    phase_sample!(happair1, happair2, compressed_Hunique, [seen])

Phases `happair1` and `happair2` window-by-window using a heuristic strategy.
`happair1[w]` can be mutated to an equivalent haplotype in window `w` such
that overall windows the number of haplotype switches is minimzed.

# Arguments
- `happair1`: `happair1[w]` is the haplotype index for strand1 in window `w`
- `happair2`: `happair2[w]` is the haplotype index for strand2 in window `w`
- `compressed_Hunique`: A `CompressedHaplotypes` object

# Optional storage argument
- `seen`: Preallocated storage container
"""
function phase_sample!(
    happair1::AbstractVector{<:Integer},
    happair2::AbstractVector{<:Integer},
    compressed_Hunique::CompressedHaplotypes,
    # preallocated items
    survivors1::AbstractVector{<:Integer}=Int32[],
    survivors2::AbstractVector{<:Integer}=Int32[],
    )

    windows = length(happair1)
    windows == length(happair2) || error("happair1 and happair2" *
        " have different lengths.")
    haplotypes = nhaplotypes(compressed_Hunique)
    lifespan1 = lifespan2 = 1 # counter to track survival time

    # get first window's optimal haplotypes
    h1 = happair1[1]
    h2 = happair2[1]
    store!(survivors1, get(compressed_Hunique.CW_typed[1].hapmap, h1, h1))
    store!(survivors2, get(compressed_Hunique.CW_typed[1].hapmap, h2, h2))

    @inbounds for w in 2:windows
        # get current window's best haplotypes
        h1 = happair1[w]
        h2 = happair2[w]
        h1set = get(compressed_Hunique.CW_typed[w].hapmap, h1, h1)
        h2set = get(compressed_Hunique.CW_typed[w].hapmap, h2, h2)

        # heuristic to decide whether cross-over is better
        # A   B      A   B
        # |   |  or    X
        # C   D      C   D
        AC = intersect_size_sorted(survivors1, h1set)
        BD = intersect_size_sorted(survivors2, h2set)
        AD = intersect_size_sorted(survivors1, h2set)
        BC = intersect_size_sorted(survivors2, h1set)
        crossed = false
        if AC + BD < AD + BC
            crossed = true
        end

        # Prune survivors. If there are none, record last survivor into
        # happair for all previous windows and reset survivors to current
        # window's haplotypes
        if crossed
            if AD == 0 # no survivors
                happair1[(w - lifespan1):(w - 1)] .= survivors1[1]
                store!(survivors1, h2set)
                lifespan1 = 1
            else
                intersect_lange!(survivors1, h2set)
                lifespan1 += 1
            end
            if BC == 0 # no survivors
                happair2[(w - lifespan2):(w - 1)] .= survivors2[1]
                store!(survivors2, h1set)
                lifespan2 = 1
            else
                intersect_lange!(survivors2, h1set)
                lifespan2 += 1
            end
        else
            if AC == 0 # no survivors
                happair1[(w - lifespan1):(w - 1)] .= survivors1[1]
                store!(survivors1, h1set)
                lifespan1 = 1
            else
                intersect_lange!(survivors1, h1set)
                lifespan1 += 1
            end
            if BD == 0 # no survivors
                happair2[(w - lifespan2):(w - 1)] .= survivors2[1]
                store!(survivors2, h2set)
                lifespan2 = 1
            else
                intersect_lange!(survivors2, h2set)
                lifespan2 += 1
            end
        end
    end

    # treat last few windows separately since intersection may not become empty
    happair1[(windows - lifespan1 + 1):windows] .= survivors1[1]
    happair2[(windows - lifespan2 + 1):windows] .= survivors2[1]

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
function intersect!(
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
function intersect!(
    v::AbstractVector{<:Integer},
    u::Integer,
    seen::AbstractSet
    )
    keep = u in v
    empty!(v)
    keep && push!(v, u)
    nothing
end

"""
    intersect_size(v::AbstractVector, u::AbstractVector, seen::BitSet=BitSet())

Computes the size of `v ∩ u` in place. Assumes `v` is usually smaller than `u`
and each element in `v` is unique.

# Arguments
- `v`: An integer vector
- `u`: An integer vector
- `seen`: Preallocated storage container
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
    for i in eachindex(v)
        @inbounds v[i] ∈ seen && (s += 1)
    end
    return s
end
intersect_size(v::AbstractVector, u::Integer, seen) = u in v

"""
    store!(v, u)

Deletes everything in `v` and saves each element of `u` to `v`.
"""
function store!(v::AbstractVector, u)
    empty!(v)
    for i in u
        push!(v, i)
    end
    return nothing
end

"""
    intersect_sorted!(v::AbstractVector, u::AbstractVector)

Computes `v ∩ u` in place and stores result in `v`. `v` and `u` is assumed sorted.
Repeated elements is allowed. 
"""
function intersect_sorted!(
    v::AbstractVector{<:Integer}, 
    u::AbstractVector{<:Integer}
    )
    lv = length(v)
    lu = length(u)
    i  = j = 1
    @inbounds while i ≤ lv && j ≤ lu
        if v[i] == u[j]
            i += 1
            j += 1
        elseif v[i] > u[j]
            j += 1
        else
            deleteat!(v, i)
            lv -= 1
        end
    end
    # handle remainders
    while i ≤ lv
        deleteat!(v, i)
        lv -= 1
    end
    return nothing
end
function intersect_sorted!(
    v::AbstractVector{<:Integer},
    u::Integer
    )
    keep = u in v
    empty!(v)
    keep && push!(v, u)
    nothing
end

"""
    intersect_size_sorted(v::AbstractVector, u::AbstractVector)

Computes the size of `v ∩ u` in place. `v` and `u` is assumed sorted.
Repeated elements is allowed. 

# Arguments
- `v`: An integer vector
- `u`: An integer vector
- `seen`: Preallocated storage container
"""
function intersect_size_sorted(
    v::AbstractVector{<:Integer}, 
    u::AbstractVector{<:Integer}
    )
    lv = length(v)
    lu = length(u)
    s = 0
    i = 1
    j = 1
    @inbounds while i ≤ lv && j ≤ lu
        if v[i] == u[j]
            s += 1
            i += 1
            j += 1
        elseif v[i] > u[j]
            j += 1
        else
            i += 1
        end
    end
    return s
end
@inline intersect_size_sorted(v::AbstractVector, u::Integer) = u in v

function intersect_lange!(v::Vector{T}, u::Vector{T}) where T <: Integer
    lv = length(v)
    lu = length(u)
    s = 0
    i = 1
    j = 1
    @inbounds while i ≤ lv && j ≤ lu
        if v[i] == u[j]
            s += 1
            v[s] = v[i]
            i += 1
            j += 1
        elseif v[i] > u[j]
            j += 1
        else
            i += 1
        end
    end
    deleteat!(v, s+1:lv)
    return nothing
end
intersect_lange!(v::AbstractVector{<:Integer}, u::Integer) = 
    intersect_sorted!(v, u)
