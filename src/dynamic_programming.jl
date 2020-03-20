"""
Helper function to calculate the difference between 2 tuples. 

# Inputs 
- `pair1`
- `pair2`

# Optional Inputs
- `λ`: Error each switch contributes. Defaults to 1.0

# Examples
- `pair_error((1, 2), (2, 3)) = pair_error((2, 1), (2, 3)) = λ`
- `pair_error((2, 5), (5, 2)  = 0` 
- `pair_error((1, 2), (3, 4)) = 2λ`
"""
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

"""
Finds the optimal sequence of haplotype pairs across all windows 
such that number of switch points is minimized. 

# Inputs
- `haplotype_set`: A person's possible haplotype pairs in each window. 

# Optional input:
- `λ`: Error each switch contributes. Defaults to 1.0

# Output:
- `sol_path`: Optimal sequence of haplotype pairs in each window
- `memory`: Vector of dictionary storing the optimal error for each haplotype pair in each window
- `path_err`: Error for each window induced by `sol_path`
- `best_err`: Osverall error induced by `sol_path`. Equals λ times number of switch points. 
"""
function connect_happairs(
    haplotype_set::Vector{Vector{T}};
    λ::Float64 = 1.0
    ) where T <: Tuple{Int, Int}

    # allocate working arrays
    windows  = length(haplotype_set)
    sol_path = Vector{T}(undef, windows)
    memory   = [Dict{T, Tuple{Float64, T}}() for i in 1:windows]

    # computational routine
    best_err = connect_happairs!(sol_path, memory, haplotype_set, λ = λ)

    return sol_path, memory, best_err
end

"""
In-place version of `connect_happairs`. 

# Inputs
- `sol_path`: Optimal sequence of haplotype pairs in each window
- `memory`: Vector of dictionary storing the optimal error for each haplotype pair in each window
- `path_err`: Error for each window induced by `sol_path`
- `haplotype_set`: A vector of vectors. `haplotype_set[1]` stores all pairs of haplotypes in window 1 in a vector, and so on. 
- `λ`: Error each switch contributes. Defaults to 1.0

# Output
- `best_err`: Overall error induced by `sol_path`. Equals λ times number of switch points. 
"""
function connect_happairs!(
    sol_path::Vector{T},
    memory::Vector{Dict{T, P}},
    haplotype_set::Vector{Vector{T}};
    λ::Float64 = 1.0,
    ) where {T <: Tuple{Int, Int}, P <: Tuple{Float64, T}}

    windows = length(haplotype_set)
    empty!.(memory) # reset storage

    # base case: last window induces no error and connects to nothing
    for pair in haplotype_set[windows]
        memory[windows][pair] = (0.0, (0, 0))
    end

    # search for best haplotype pair in each window bottom-up 
    for w in Iterators.reverse(1:(windows - 1)), happair in haplotype_set[w]
        # search all pairs in next window
        best_err = Inf
        best_next_pair = (0, 0)
        for pair in haplotype_set[w + 1]
            err = pair_error(happair, pair) + memory[w + 1][pair][1]
            if err < best_err
                best_err = err
                best_next_pair = pair
            end
        end
        memory[w][happair] = (best_err, best_next_pair)
    end

    # find best starting point
    best_err   = Inf
    best_start = (0, 0) 
    for (key, val) in memory[1]
        if val[1] < best_err
            best_start = key
            best_err   = val[1]
        end
    end

    # find best solution path by tracing `memory`
    sol_path[1] = best_start
    for w in 2:windows
        prev_pair   = sol_path[w - 1]
        sol_path[w] = memory[w - 1][prev_pair][2]
    end

    return best_err
end
