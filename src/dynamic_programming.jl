"""
By exchanging entries of 2 vectors at the same position, this 
function minimize the number of breakpoints in both strands, 
and ruturn that number.

# Arguments:
- `n`: Starting position. Flipping will always happen at `n - 1`
- `seq1`: First sequence
- `seq2`: Second sequence
- `flip_idx`: A preallocated vector that will store the index where flipping occurred
- `intermediates`: vector storing intermediate results (for memoization). `-1` indicates not solved.

# Example 1
x = [0 1 1 1]
y = [1 0 0 0]
We can exchange position 1:
x = [1 1 1 1]
y = [0 0 0 0]
then breakpoint = 0

# Example 2
x = [0 1 0 0 1]
y = [1 0 1 0 0]
we can exchange position 2:
x = [0 0 0 0 1]
y = [1 1 1 0 0]
then breakpoint = 2 (1 at end of x and another at middle of y)
"""
function binary_flip!(n, seq1, seq2, flip_idx,
    intermediates::Vector{Int} = [-1 for i in 1:(length(seq1) + 1)]
    )
    # quick lookup
    if intermediates[n] != -1
        return intermediates[n]
    end

    if n > length(seq1)
        return 0
    elseif n == 1
        # only flip bits in previous position
        return binary_flip!(n + 1, seq1, seq2, flip_idx, intermediates)
    elseif seq1[n] == seq2[n]
        # don't flip and only calculate error
        return binary_flip!(n + 1, seq1, seq2, flip_idx, intermediates) + !(seq1[n] == seq1[n - 1]) + !(seq2[n] == seq2[n - 1])
    else
        # calculate error of flip/noflip by recursion
        yesflip = binary_flip!(n + 1, seq1, seq2, flip_idx, intermediates) + !(seq1[n] == seq2[n - 1]) + !(seq2[n] == seq1[n - 1])
        noflip  = binary_flip!(n + 1, seq1, seq2, flip_idx, intermediates) + !(seq1[n] == seq1[n - 1]) + !(seq2[n] == seq2[n - 1])
        # println("n = $n, yesflip = $yesflip, noflip = $noflip, min = $(min(yesflip, noflip))")
        # store intermediate results for later retrival
        intermediates[n] = min(yesflip, noflip)
        if yesflip < noflip
            # record flipping location, flip 2 sequence at previous location, 
            flip_idx[n - 1] = true
            seq1[n - 1], seq2[n - 1] = seq2[n - 1], seq1[n - 1]
            # println("flipped $(n - 1)!")
        end
        return min(yesflip, noflip)
    end
end
binary_flip!(seq1, seq2, flip) = binary_flip!(2, seq1, seq2, flip) #start at position 2

function has_intersect!(c::BitVector, a::BitVector, b::BitVector)
    c .= a .& b
    return any(c)
end

"""
Essentially the same code as `binary_flip!`. The difference is, instead 
of exchanging 0s and 1s, we exchange sets.  
"""
function set_flip!(
    n::Int, 
    vector_set1::Vector{BitVector}, 
    vector_set2::Vector{BitVector}, 
    flip_idx::BitVector,
    intermediates::Vector{Int} = [-1 for i in 1:(length(vector_set1) + 1)],
    storage::BitVector = falses(length(vector_set1[1]))
    )
    # quick lookup
    if intermediates[n] != -1
        return intermediates[n]
    end

    if n > length(vector_set1)
        return 0
    elseif n == 1
        # only flip bits in previous position
        return set_flip!(n + 1, vector_set1, vector_set2, flip_idx, intermediates)
    elseif vector_set1[n] == vector_set2[n]
        # don't flip and only calculate error
        return set_flip!(n + 1, vector_set1, vector_set2, flip_idx, intermediates) + 
                !has_intersect!(storage, vector_set1[n], vector_set1[n - 1]) +
                !has_intersect!(storage, vector_set2[n], vector_set2[n - 1])
    else
        # calculate error of flip/noflip by recursion
        yesflip = set_flip!(n + 1, vector_set1, vector_set2, flip_idx, intermediates) + 
                    !has_intersect!(storage, vector_set1[n], vector_set2[n - 1]) + 
                    !has_intersect!(storage, vector_set2[n], vector_set1[n - 1])
        noflip  = set_flip!(n + 1, vector_set1, vector_set2, flip_idx, intermediates) + 
                    !has_intersect!(storage, vector_set1[n], vector_set1[n - 1]) + 
                    !has_intersect!(storage, vector_set2[n], vector_set2[n - 1])

        # store intermediate results for later retrival
        intermediates[n] = min(yesflip, noflip)
        if yesflip < noflip
            # record flipping location, flip 2 sequence at previous location, 
            flip_idx[n - 1] = true
            vector_set1[n - 1], vector_set2[n - 1] = vector_set2[n - 1], vector_set1[n - 1]
        end
        return min(yesflip, noflip)
    end
end
set_flip!(vector_set1, vector_set2, flip) = set_flip!(2, vector_set1, vector_set2, flip) #start at position 2

"""
Finds the optimal sequence of haplotype pairs across all windows 
such that number of switch points is minimized. 

# Inputs
- `w`: Current window
- `happair`: Haplotype pair in current window being considered. 
- `haplotype_set`: A vector of vectors. `haplotype_set[1]` stores all pairs of haplotypes in window 1 in a vector, and so on. 
- `λ`: Error associated with 1 mismatch when comparing 2 pairs of haplotypes. e.g. (h1, h2) vs (h1, h3) have error λ. 
"""
function connect_happairs(haplotype_set::Vector{Vector{T}}; λ::Float64 = 1.0) where T <: Tuple{Int, Int}
    windows  = length(haplotype_set)
    memory   = [Dict{T, Float64}() for i in 1:windows]
    sol_path = Vector{T}(undef, windows)

    best_err  = Inf
    best_pair1 = (0, 0)
    best_pair2 = (0, 0)
    for happair in haplotype_set[1], pair in haplotype_set[2]
        err = pair_error(happair, pair) + connect_happairs(2, pair, haplotype_set, λ = λ, memory = memory, solution_path = sol_path)
        if err < best_err
            best_pair1 = happair
            best_pair2 = pair
            best_err   = err
        end
    end

    # save best pair in first 2 windows
    sol_path[1] = best_pair1
    sol_path[2] = best_pair2

    return sol_path, memory, best_err
end

"""
# Inputs
- `w`: Current window
- `happair`: Haplotype pair in current window being considered. 
- `haplotype_set`: A vector of vectors. `haplotype_set[1]` stores all pairs of haplotypes in window 1 in a vector, and so on. 
- `λ`: Error associated with 1 mismatch when comparing 2 pairs of haplotypes. e.g. (h1, h2) vs (h1, h3) have error λ. 
"""
function connect_happairs(
    w::Int,
    happair::T, 
    haplotype_set::Vector{Vector{T}};
    λ::Float64 = 1.0,
    memory = [Dict{T, Float64}() for i in 1:length(haplotype_set)],
    solution_path = Vector{T}(undef, length(haplotype_set))
    ) where T <: Tuple{Int, Int}

    if haskey(memory[w], happair)
        return memory[w][happair] # quick lookup
    elseif w == length(haplotype_set)
        return 0 # last window contributes no extra error
    else
        # recursion: solve next subtree 
        best_err = Inf
        best_next_pair = (0, 0)
        for pair in haplotype_set[w + 1]
            err = pair_error(happair, pair) + connect_happairs(w + 1, pair, haplotype_set, λ = λ, memory = memory, solution_path = solution_path)
            if err < best_err
                best_next_pair = pair
                best_err = err
            end
        end

        # record best error and solution path
        memory[w][happair] = best_err
        solution_path[w + 1] = best_next_pair
        return best_err
    end
end

function pair_error(pair1::T, pair2::T; λ::Real = 1.0) where T <: Tuple{Int, Int}
    difference = zero(eltype(λ))
    if pair1[1] != pair2[1] && pair1[1] != pair2[2]
        difference += one(eltype(λ))
    end
    if pair1[2] != pair2[1] && pair1[2] != pair2[2]
        difference += one(eltype(λ))
    end
    return λ * difference
end
