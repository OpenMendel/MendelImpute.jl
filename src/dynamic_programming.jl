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

