"""
Helper function to calculate the squared Hamming difference between 2 unordered pair of integers.

# Inputs 
- `pair1`: tuple (a, b) where a, b are integers
- `pair2`: tuple (c, d) where c, d are integers

# Optional Inputs
- `λ`: Error each switch contributes. Defaults to 1.0

# Examples
- `pair_error((1, 2), (2, 3)) = pair_error((2, 1), (2, 3)) = λ`
- `pair_error((2, 5), (5, 2)  = 0` 
- `pair_error((1, 2), (3, 4)) = 4λ`
"""
function pair_error(pair1::T, pair2::T; λ::Real = 1.0) where T <: Tuple{Int, Int}
    # parallel connections
    # a b
    # | |
    # c d
    parallel_diff = (pair1[1] != pair2[1]) + (pair1[2] != pair2[2])
    # a b
    #  X
    # c d
    crossover_diff = (pair1[1] != pair2[2]) + (pair1[2] != pair2[1])
    return λ * abs2(min(parallel_diff, crossover_diff))
end

"""
Finds the optimal sequence of haplotype pairs across all windows 
such that number of switch points is minimized. Windows with
too little typed SNPs will inherit haplotype pairs from the closest
window containing feasible haplotype pairs. 

# Inputs
- `haplotype_set`: A person's possible haplotype pairs in each window. `(-1, -1)` indicates too little typed snps. 

# Optional input:
- `λ`: Error each switch contributes. Defaults to 1.0

# Output:
- `sol_path`: Optimal sequence of haplotype pairs indices (in the complete haplotype pool) in each window
- `memory`: Vector of dictionary storing the optimal error for each haplotype pair in each window
- `path_err`: Error for each window induced by `sol_path`
- `best_err`: Osverall error induced by `sol_path`. Equals λ times number of switch points. 
"""
function connect_happairs!(
    haplotype_set::Vector{Vector{T}};
    λ::Float64 = 1.0
    ) where T <: Tuple{Int, Int}

    # allocate working arrays
    windows  = length(haplotype_set)
    sol_path = Vector{T}(undef, windows)
    next_pair = [Int[] for i in 1:windows]
    subtree_err = [Float64[] for i in 1:windows]

    # computational routine
    best_err = connect_happairs!(sol_path, next_pair, subtree_err, haplotype_set, λ = λ)

    return sol_path, next_pair, subtree_err, best_err
end

"""
Finds the optimal sequence of haplotype pairs across all windows 
such that number of switch points is minimized. Windows with
too little typed SNPs will inherit haplotype pairs from the closest
window containing feasible haplotype pairs. 

# Inputs
- `sol_path`: Optimal sequence of haplotype pairs indices (in the complete haplotype pool) in each window
- `memory`: Vector of dictionary storing the optimal error for each haplotype pair in each window
- `path_err`: Error for each window induced by `sol_path`
- `haplotype_set`: A vector of vectors. `haplotype_set[1]` stores all pairs of haplotypes in window 1 in a vector, and so on. 
- `λ`: Error each switch contributes. Defaults to 1.0

# Optional input:
- `λ`: Error each switch contributes. Defaults to 1.0

# Output
- `best_err`: Overall error induced by `sol_path`. Equals λ times number of switch points. 
"""
function connect_happairs!(
    sol_path::Vector{T},
    next_pair::Vector{Vector{Int}}, 
    subtree_err::Vector{Vector{Float64}},
    haplotype_set::Vector{Vector{T}};
    λ::Float64 = 1.0,
    ) where T <: Tuple{Int, Int}

    windows = length(haplotype_set)

    # reset storage
    empty!.(next_pair) 
    empty!.(subtree_err)

    # base case: last window induces no error and connects to nothing
    @inbounds for pair in haplotype_set[windows]
        push!(next_pair[windows], 0)
        push!(subtree_err[windows], 0.0)
    end

    # search for best haplotype pair in each window bottom-up 
    @inbounds for w in Iterators.reverse(1:(windows - 1))
        win_best_err = Inf

        # first pass to compute each pair's optimal error
        for (j, happair) in enumerate(haplotype_set[w])
            # search all pairs in next window
            best_err = Inf
            best_next_pair = 0
            for (i, nextpair) in enumerate(haplotype_set[w + 1])
                err = pair_error(happair, nextpair) + subtree_err[w + 1][i]
                if err < best_err
                    best_err = err
                    best_next_pair = i
                end
            end
            if best_err < win_best_err
                win_best_err = best_err
            end
            push!(subtree_err[w], best_err)
            push!(next_pair[w], best_next_pair)
        end

        # remove pairs that are much worse than best error, reducing dp search space
        tol = win_best_err # remove all suboptimal pairs
        for (i, err) in enumerate(subtree_err[w])
            if err > tol
                deleteat!(subtree_err[w], i)
                deleteat!(next_pair[w], i)
                deleteat!(haplotype_set[w], i)
            end
        end
    end

    # find best solution path by forward-tracing
    best_err, cur_idx = findmin(subtree_err[1])
    @inbounds for w in 1:windows
        sol_path[w] = haplotype_set[w][cur_idx]
        cur_idx = next_pair[w][cur_idx]
    end

    return best_err
end

# function nearest_window_with_sufficient_typed_snps(
#     cur_win::Int, 
#     haplotype_set::Vector{Vector{T}},
#     originally_empty::BitVector
#     ) where T <: Tuple{Int, Int}

#     # quick return
#     if haplotype_set[cur_win][1] != (-1, -1)
#         return cur_win
#     end

#     total_windows = length(haplotype_set)
#     right_windows = total_windows - cur_win
#     left_windows  = cur_win - 1

#     #extend to the right until we find a window with feasible happair
#     right = typemax(Int)
#     for w in 1:right_windows
#         next_win = cur_win + w
#         if haplotype_set[next_win][1] != (-1, -1) && !originally_empty[next_win]
#             right = next_win
#             break
#         end
#     end

#     #extend to the left until we find a window with feasible happair
#     left = typemax(Int)
#     for w in 1:left_windows
#         next_win = cur_win - w
#         if haplotype_set[next_win][1] != (-1, -1) && !originally_empty[next_win]
#             left = next_win
#             break
#         end
#     end 

#     right == left == typemax(Int) && error("No window contains feasible haplotype pairs!")
    
#     if abs(cur_win - right) < abs(cur_win - left)
#         return right
#     else
#         return left
#     end
# end
