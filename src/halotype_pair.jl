"""
    unique_haplotypes(H, width, trans)

For each window, finds unique haplotype indices stored in the columns of H and 
saves a mapping vector of unique columns of H. See `UniqueHaplotypeMaps` data 
structure for examples. 

# Input
* `H`: An `p x d` reference panel of haplotypes within a genomic window. 
* `width`: The window width 
* `trans`: Orientation of `H`. 'T' means columns of `H` are a haplotype vectors. 'N' means rows of `H` are. 

# Output
* `hapset`: Data structure for keeping track of unique haplotypes in each window. 
"""
function unique_haplotypes(
    H::AbstractMatrix,
    width::Int,
    trans::Char;
    flankwidth::Int = round(Int, 0.1width)
    )

    if trans == 'N'
        dim = 1
    elseif trans == 'T'
        dim = 2
    else
        error("trans can only be 'N' or 'T' but was $dim" )
    end

    p, d    = size(H)
    windows = floor(Int, p / width)
    hapset  = UniqueHaplotypeMaps(windows, d)

    # record unique haplotypes and mappings window by window (using flanking windows)
    # first  1/3: ((w - 2) * width + 1):((w - 1) * width)
    # middle 1/3: ((w - 1) * width + 1):(      w * width)
    # last   1/3: (      w * width + 1):((w + 1) * width)
    Threads.@threads for w in 1:windows
        if w == 1
            cur_range = 1:(width + flankwidth)
        elseif w == windows
            cur_range = ((windows - 1) * width - flankwidth + 1):p
        else
            cur_range = ((w - 1) * width - flankwidth + 1):(w * width + flankwidth)
        end

        H_cur_window = view(H, cur_range, :)
        hapset.hapmap[w] = groupslices(H_cur_window, dim)
        hapset.uniqueindex[w] = unique(hapset.hapmap[w])
        hapset.range[w] = cur_range
    end

    return hapset
end

"""
    compute_optimal_halotype_set(X, H, width, verbose)

Computes the optimal haplotype pair for each person in each window, then computes
a set of haplotypes that matches the optimal haplotype pair in the current window.

# Input 
+ `X`: Target genotype matrix with missing entries. Each column is a person's genotype
+ `H`: Reference haplotype panels, each column is a haplotype. 
+ `width`: The width of each window
+ `verbose`: boolean indicating whether to print intermediate results. 

# Output
+ `optimal_haplotypes`: where optimal_haplotypes[i] is a `OptimalHaplotypeSet` recording all
redundant haplotypes that matches the optimal haplotypes in each window for person i. 
"""
function compute_optimal_halotype_set(
    X::AbstractMatrix,
    H::AbstractMatrix;
    width::Int = 400,
    flankwidth::Int = round(Int, 0.1width),
    verbose::Bool = true,
    fast_method::Bool = false,
    )

    # define some constants
    snps, people = size(X)
    haplotypes = size(H, 2)
    windows = floor(Int, snps / width)
    threads = Threads.nthreads()
    T = Float32

    # Each person stores a vector of redundant haplotypes matching the optimal one for each window
    if fast_method
        redundant_haplotypes = [OptimalHaplotypeSet(windows, haplotypes) for i in 1:people]
    else
        redundant_haplotypes = [[Tuple{Int, Int}[] for i in 1:windows] for j in 1:people]
    end

    # get unique haplotype indices and maps for each window
    Hunique = unique_haplotypes(H, width, 'T', flankwidth = flankwidth)

    # In first window, calculate optimal haplotype pair among unique haplotypes
    pmeter   = Progress(windows, 5, "Computing optimal haplotype pairs...")
    happairs = [Tuple{Int, Int}[] for i in 1:people] # tracks unique haplotype pairs in a window
    hapscore = zeros(T, people)
    num_uniq = length(Hunique.uniqueindex[1])
    M        = zeros(T, num_uniq, num_uniq)
    N        = zeros(T, people, num_uniq)
    cur_range = Hunique.range[1]
    Hwork_tmp = convert(Matrix{T}, @view(H[cur_range, Hunique.uniqueindex[1]]))
    Xwork_tmp = @view(X[cur_range, :])
    haploimpute!(Xwork_tmp, Hwork_tmp, M, N, happairs, hapscore)
    compute_redundant_haplotypes!(redundant_haplotypes, Hunique, happairs, 1, fast_method=fast_method)
    next!(pmeter)

    # allocate working arrays for remaining windows (beware: window 1's cur_range length is different than these)
    M = Vector{Matrix{T}}(undef, threads)
    N = Vector{ElasticArray{T}}(undef, threads)
    Xwork       = Vector{AbstractMatrix{Union{eltype(X), Missing}}}(undef, threads)
    Xwork_float = Vector{Matrix{T}}(undef, threads)
    Hwork_float = Vector{ElasticArray{T}}(undef, threads)
    for id in 1:threads 
        win = id + 1 #avoid window 1
        cur_range = Hunique.range[win]
        num_uniq  = length(cur_range)
        Hwork_float[id] = ElasticArray{T}(@view(H[cur_range, Hunique.uniqueindex[win]]))
        Xwork_float[id] = zeros(T, num_uniq, people)
        M[id] = zeros(T, num_uniq, num_uniq)
        N[id] = ElasticArray{T}(undef, people, num_uniq)
    end
    happairs = [[Tuple{Int, Int}[] for i in 1:people] for _ in 1:threads]
    hapscore = [zeros(T, people) for _ in 1:threads]

    # loop through each window
    Threads.@threads for w in 2:(windows - 1)
        id = Threads.threadid()

        # sync working arrays with current window's data
        next_dim = length(Hunique.uniqueindex[w])
        Xwork[id] = view(X, Hunique.range[w], :)
        size(M[id], 1) != next_dim && (M[id] = zeros(T, next_dim, next_dim)) # Julia can't resize matrix, at least not till 2.0
        resize_and_sync!(Hwork_float[id], N[id], Hunique.uniqueindex[w], Hunique.range[w], H)

        # Calculate optimal haplotype pair among unique haplotypes
        haploimpute!(Xwork[id], Hwork_float[id], M[id], N[id], happairs[id], hapscore[id], Xfloat=Xwork_float[id])

        # find all haplotypes matching the optimal haplotype pairs
        compute_redundant_haplotypes!(redundant_haplotypes, Hunique, happairs[id], w, fast_method=fast_method)

        # update progress
        next!(pmeter)
    end

    # last window reallocate everything 
    last_range  = Hunique.range[end]
    num_uniq    = length(Hunique.uniqueindex[end])
    Hwork_float = convert(Matrix{T}, @view(H[last_range, Hunique.uniqueindex[end]]))
    Xwork       = @view(X[last_range, :])
    M           = zeros(T, num_uniq, num_uniq)
    N           = zeros(T, people, num_uniq)
    haploimpute!(Xwork, Hwork_float, M, N, happairs[1], hapscore[1])
    compute_redundant_haplotypes!(redundant_haplotypes, Hunique, happairs[1], windows, fast_method=fast_method)
    next!(pmeter)

    return redundant_haplotypes
end

function compute_optimal_halotype_pair(
    X::AbstractMatrix{Union{Missing, T}},
    H::AbstractMatrix;
    width::Int = 400,
    flankwidth::Int = round(Int, 0.1width),
    fast_method::Bool = false,
    ) where T <: Real

    # define some constants
    snps, people = size(X)
    haplotypes = size(H, 2)
    windows = floor(Int, snps / width)
    threads = Threads.nthreads()

    # get unique haplotype indices and maps for each window
    Hunique = unique_haplotypes(H, width, 'T', flankwidth = flankwidth)

    # stores optimal haplotype pairs in each window
    optimal_happairs = [[Tuple{Int, Int}[] for i in 1:windows] for j in 1:people]

    # allocate working arrays
    happairs    = [Tuple{Int, Int}[] for i in 1:people] # tracks unique haplotype pairs in a window
    hapscore    = zeros(T, people)
    num_uniq    = length(Hunique.uniqueindex[1])
    M           = zeros(T, num_uniq, num_uniq)
    N           = ElasticArray{T}(undef, people, num_uniq) # array type that allows rescaling last dim 
    pmeter      = Progress(windows, 5, "Computing optimal haplotype pairs...")

    # In first window, calculate optimal haplotype pair among unique haplotypes
    cur_range   = Hunique.range[1]
    Hwork_tmp   = convert(Matrix{T}, @view(H[cur_range, Hunique.uniqueindex[1]]))
    Xwork_tmp   = @view(X[cur_range, :])
    haploimpute!(Xwork_tmp, Hwork_tmp, M, N, happairs, hapscore)
    for k in 1:people, happair in happairs[k]
        Hi_uniqueidx, Hj_uniqueidx = happair
        Hi_idx = Hunique.uniqueindex[1][Hi_uniqueidx]
        Hj_idx = Hunique.uniqueindex[1][Hj_uniqueidx]
        push!(optimal_happairs[k][1], (Hi_idx, Hj_idx))
    end
    next!(pmeter)

    # new resizable working arrays for remaining windows since window 1's size may be different
    M = Vector{Matrix{T}}(undef, threads)
    N = Vector{ElasticArray{T}}(undef, threads)
    Xwork       = Vector{AbstractMatrix{Union{T, Missing}}}(undef, threads)
    Xwork_float = Vector{Matrix{T}}(undef, threads)
    Hwork_float = Vector{ElasticArray{T}}(undef, threads)
    for id in 1:threads 
        win = id + 1 #avoid window 1
        cur_range = Hunique.range[win]
        num_uniq  = length(cur_range)
        Hwork_float[id] = ElasticArray{T}(@view(H[cur_range, Hunique.uniqueindex[win]]))
        Xwork_float[id] = zeros(T, num_uniq, people)
        M[id] = zeros(T, num_uniq, num_uniq)
        N[id] = ElasticArray{T}(undef, people, num_uniq)
    end
    happairs = [[Tuple{Int, Int}[] for i in 1:people] for _ in 1:threads]
    hapscore = [zeros(T, people) for _ in 1:threads]

    # loop through each window
    Threads.@threads for w in 2:(windows - 1)
        id = Threads.threadid()

        # sync working arrays with current window's data
        next_dim = length(Hunique.uniqueindex[w])
        Xwork[id] = view(X, Hunique.range[w], :)
        size(M[id], 1) != next_dim && (M[id] = zeros(T, next_dim, next_dim)) # Julia can't resize matrix, at least not till 2.0
        resize_and_sync!(Hwork_float[id], N[id], Hunique.uniqueindex[w], Hunique.range[w], H)

        # Calculate optimal haplotype pair among unique haplotypes
        haploimpute!(Xwork[id], Hwork_float[id], M[id], N[id], happairs[id], hapscore[id], Xfloat=Xwork_float[id])

        # store current window's optimal happairs
        for k in 1:people, happair in happairs[id][k]
            Hi_uniqueidx, Hj_uniqueidx = happair
            Hi_idx = Hunique.uniqueindex[w][Hi_uniqueidx]
            Hj_idx = Hunique.uniqueindex[w][Hj_uniqueidx]
            push!(optimal_happairs[k][w], (Hi_idx, Hj_idx))
        end

        # update progress
        next!(pmeter)
    end

    # last window reallocate everything 
    last_range  = Hunique.range[end]
    num_uniq    = length(Hunique.uniqueindex[end])
    Hwork_float = convert(Matrix{T}, @view(H[last_range, Hunique.uniqueindex[end]]))
    Xwork       = @view(X[last_range, :])
    M           = zeros(T, num_uniq, num_uniq)
    N           = zeros(T, people, num_uniq)
    haploimpute!(Xwork, Hwork_float, M, N, happairs[1], hapscore[1])
    for k in 1:people, happair in happairs[1][k]
        Hi_uniqueidx, Hj_uniqueidx = happair
        Hi_idx = Hunique.uniqueindex[windows][Hi_uniqueidx]
        Hj_idx = Hunique.uniqueindex[windows][Hj_uniqueidx]
        push!(optimal_happairs[k][windows], (Hi_idx, Hj_idx))
    end
    next!(pmeter)

    return optimal_happairs
end

# function compute_optimal_halotype_set_prephased(
#     X::AbstractMatrix{Union{Missing, T}},
#     H::AbstractMatrix;
#     width::Int    = 400,
#     flankwidth::Int = round(Int, 0.1width),
#     ) where T <: Real
    
#     # declare some constants
#     snps = size(X, 1)
#     people = Int(size(X, 2) / 2)
#     haplotypes = size(H, 2)
#     windows = floor(Int, snps / width)

#     # get unique haplotype indices and maps for current window
#     Hunique = unique_haplotypes(H, width, 'T', flankwidth = flankwidth)
#     hapset  = [OptimalHaplotypeSet(windows, haplotypes) for i in 1:people]

#     for w in 1:windows, i in 1:people
#         cur_range = Hunique.range[w]
#         Hi_unique = Hunique.uniqueindex[w]
#         best_err1 = typemax(eltype(H))
#         best_err2 = typemax(eltype(H))

#         # loop through 2 genotype strands to find best matching haplotype
#         Xi1 = @view(X[cur_range, 2i - 1])
#         Xi2 = @view(X[cur_range, 2i])
#         Hi1 = 0
#         Hi2 = 0
#         for j in Hi_unique
#             Hi = @view(H[cur_range, j])
#             # strand1
#             err = euclidean_skipmissing(Xi1, Hi)
#             if err < best_err1
#                 Hi1 = j
#                 best_err1 = err
#             end
#             # strand2
#             err = euclidean_skipmissing(Xi2, Hi)
#             if err < best_err2
#                 Hi2 = j
#                 best_err2 = err
#             end
#         end
        
#         # find all haplotypes that matches the unique one 
#         mapping = Hunique.hapmap[w]
#         for j in 1:haplotypes
#             mapping[j] == Hi1 && (hapset[i].strand1[w][j] = true)
#             mapping[j] == Hi2 && (hapset[i].strand2[w][j] = true)
#         end
#     end

#     return hapset
# end

# TODO: possible type instability
# function euclidean_skipmissing(
#     x::AbstractVector{Union{T, Missing}}, 
#     y::AbstractVector{T}
#     ) where {T <: Real}
#     s = zero(T)
#     @inbounds @simd for i in eachindex(x)
#         if x[i] !== missing
#             s += (x[i] - y[i])^2
#         end
#     end
#     return s
# end

"""
Records optimal-redundant haplotypes for each window. 
"""
function compute_redundant_haplotypes!(
    redundant_haplotypes::Union{Vector{Vector{Vector{T}}}, Vector{OptimalHaplotypeSet}}, 
    Hunique::UniqueHaplotypeMaps, 
    happairs::Vector{Vector{T}}, 
    window::Int;
    fast_method::Bool = false,
    ) where T <: Tuple{Int, Int}
    
    people = length(redundant_haplotypes)
    if fast_method
        @inbounds for k in 1:people, happair in happairs[k]
            Hi_uniqueidx = happair[1]
            Hj_uniqueidx = happair[2]
            # println("person $k's optimal haplotype pairs are: $((Hi_uniqueidx, Hj_uniqueidx))")

            Hi_idx = Hunique.uniqueindex[window][Hi_uniqueidx]
            Hj_idx = Hunique.uniqueindex[window][Hj_uniqueidx]
            # println("person $k's optimal haplotype pairs are located at columns $Hi_idx and $Hj_idx in current window of H")

            # loop through all haplotypes and find ones that match either of the optimal haplotypes 
            mapping = Hunique.hapmap[window]
            redunhaps_bitvec1 = redundant_haplotypes[k].strand1[window]
            redunhaps_bitvec2 = redundant_haplotypes[k].strand2[window]
            for jj in 1:length(Hunique.hapmap[1])
                mapping[jj] == Hi_idx && (redunhaps_bitvec1[jj] = true)
                mapping[jj] == Hj_idx && (redunhaps_bitvec2[jj] = true)
            end
        end
    else
        h1_set, h2_set = Int[], Int[]
        @inbounds for k in 1:people
            for happair in happairs[k]
                Hi_uniqueidx = happair[1]
                Hj_uniqueidx = happair[2]
                # println("person $k's optimal haplotype pairs are: $((Hi_uniqueidx, Hj_uniqueidx))")

                Hi_idx = Hunique.uniqueindex[window][Hi_uniqueidx]
                Hj_idx = Hunique.uniqueindex[window][Hj_uniqueidx]
                # println("person $k's optimal haplotype pairs are located at columns $Hi_idx and $Hj_idx in current window of H")

                # loop through all haplotypes and find ones that match either of the optimal haplotypes 
                empty!(h1_set)
                empty!(h2_set)
                for (idx, hap) in enumerate(Hunique.hapmap[window])
                    hap == Hi_idx && push!(h1_set, idx)
                    hap == Hj_idx && push!(h2_set, idx)
                end

                # push all possible happair into `redundant_haplotypes` 
                for h1 in h1_set, h2 in h2_set
                    push!(redundant_haplotypes[k][window], (h1, h2))
                end
            end
            # reduce search space for dynamic programming later
            if length(redundant_haplotypes[k][window]) > 1000
                # shuffle!(redundant_haplotypes[k][window]) # THIS IS NOT THREAD SAFE!
                resize!(redundant_haplotypes[k][window], 1000)
            end
        end
    end

    return nothing
end

"""
    resize_and_sync!(Hwork, next_d, window, H, N)

Up/downsizes the dimension of `Hwork` and `N`, and copies relevant information into `Hwork`. 

# Inputs
* `Hwork`: Haplotype matrix in the current window containing only unique haplotypes. Must add/subtract columns. 
* `N`: Matrix used in the computational routine. Must add/subtract columns. 
* `next_d`: Number of unique haplotypes in next window. 
* `window`: Indices of current window. 
* `H`: Full haplotype reference panel. Each column is a haplotype
"""
function resize_and_sync!(
    Hwork::ElasticArray,
    N::ElasticArray,
    Hnext::Vector{Int},
    window::UnitRange{Int},
    H::AbstractMatrix,
    )
    next_d = length(Hnext)
    if size(N, 2) != next_d
        resize!(N, size(N, 1), next_d)
    end
    if size(Hwork, 2) != next_d
        resize!(Hwork, size(Hwork, 1), next_d)
    end
    copyto!(Hwork, view(H, window, Hnext))
end

"""
    haplopair(X, H)

Calculate the best pair of haplotypes in `H` for each individual in `X`. Assumes `X` 
does not have missing data. 

# Input
* `X`: `p x n` genotype matrix. Each column is an individual.
* `H`: `p * d` haplotype matrix. Each column is a haplotype.

# Output
* `happair`: optimal haplotype pairs. `X[:, k] ≈ H[:, happair[1][k]] + H[:, happair[2][k]]`.
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
"""
function haplopair(
    X::AbstractMatrix,
    H::AbstractMatrix
    )

    p, n     = size(X)
    d        = size(H, 2)
    M        = zeros(eltype(H), d, d)
    N        = zeros(promote_type(eltype(H), eltype(X)), n, d)
    happairs = [Tuple{Int, Int}[] for i in 1:n]
    hapscore = zeros(eltype(N), n)
    haplopair!(X, H, M, N, happairs, hapscore)

    return happairs, hapscore
end

"""
    haplopair!(X, H, M, N, happair, hapscore)

Calculate the best pair of haplotypes in `H` for each individual in `X`. Overwite
`M` by `M[i, j] = 2dot(H[:, i], H[:, j]) + sumabs2(H[:, i]) + sumabs2(H[:, j])`,
`N` by `2X'H`, `happair` by optimal haplotype pair, and `hapscore` by
objective value from the optimal haplotype pair.

# Input
* `X`: `p x n` genotype matrix. Each column is an individual.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `M`: overwritten by `M[i, j] = 2dot(H[:, i], H[:, j]) + sumabs2(H[:, i]) +
    sumabs2(H[:, j])`.
* `N`: overwritten by `n x d` matrix `2X'H`.
* `happair`: optimal haplotype pair. `X[:, k] ≈ H[:, happair[k, 1]] + H[:, happair[k, 2]]`.
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
"""
function haplopair!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    M::AbstractMatrix,
    N::AbstractMatrix,
    happairs::Vector{Vector{Tuple{Int, Int}}},
    hapscore::AbstractVector
    )

    p, n, d = size(X, 1), size(X, 2), size(H, 2)

    # assemble M (upper triangular only)
    mul!(M, Transpose(H), H)
    for j in 1:d, i in 1:(j - 1) # off-diagonal
        M[i, j] = 2M[i, j] + M[i, i] + M[j, j]
    end
    for j in 1:d # diagonal
        M[j, j] *= 4
    end

    # assemble N
    mul!(N, Transpose(X), H)
    @simd for I in eachindex(N)
        N[I] *= 2
    end

    # computational routine
    haplopair!(happairs, hapscore, M, N)

    # supplement the constant terms in objective
    @inbounds for j in 1:n
        @simd for i in 1:p
            hapscore[j] += abs2(X[i, j])
        end
    end

    return nothing
end

"""
    haplopair!(happair, hapscore, M, N)

Calculate the best pair of haplotypes pairs in the filtered haplotype panel
for each individual in `X` using sufficient statistics `M` and `N`. 

# Note
The best haplotype pairs are column indices of the filtered haplotype 
panel, and must be converted back to the indices of the actual haplotype panel
using `UniqueHaplotypeMaps`. 

# Input
* `happair`: optimal haplotype pair for each individual.
* `hapmin`: minimum offered by the optimal haplotype pair.
* `M`: `d x d` matrix with entries `M[i, j] = 2dot(H[:, i], H[:, j]) +
    sumabs2(H[:, i]) + sumabs2(H[:, j])`, where `H` is the haplotype matrix
    with haplotypes in columns. Only the upper triangular part of `M` is used.
* `N`: `n x d` matrix `2X'H`, where `X` is the genotype matrix with individuals
    in columns.
"""
function haplopair!(
    happairs::Vector{Vector{Tuple{Int, Int}}},
    hapmin::Vector{T},
    M::AbstractMatrix{T},
    N::AbstractMatrix{T},
    ) where T <: Real

    n, d = size(N)
    tol = convert(T, 3)
    fill!(hapmin, typemax(T))
    empty!.(happairs)

    @inbounds for k in 1:d, j in 1:k
        # loop over individuals
        @simd for i in 1:n
            score = M[j, k] - N[i, j] - N[i, k]

            # keep best happair (original code)
            # if score < hapmin[i]
            #     empty!(happairs[i])
            #     push!(happairs[i], (j, k))
            #     hapmin[i] = score
            # end

            # keep all happairs that are equally good
            # if score < hapmin[i]
            #     empty!(happairs[i])
            #     push!(happairs[i], (j, k))
            #     hapmin[i] = score
            # elseif score == hapmin[i]
            #     push!(happairs[i], (j, k))
            # end

            # keep happairs that within some range of best pair (but finding all of them requires a 2nd pass)
            if score < hapmin[i]
                empty!(happairs[i])
                push!(happairs[i], (j, k))
                hapmin[i] = score
            elseif score <= hapmin[i] + tol
                push!(happairs[i], (j, k))
            end

            # keep top 10 haplotype pairs
            # if score < hapmin[i]
            #     length(happairs[i]) == 10 && popfirst!(happairs[i])
            #     push!(happairs[i], (j, k))
            #     hapmin[i] = score
            # elseif score <= hapmin[i] + interval
            #     length(happairs[i]) == 10 && popfirst!(happairs[i])
            #     push!(happairs[i], (j, k))
            # end

            # keep all previous best pairs
            # if score < hapmin[i]
            #     push!(happairs[i], (j, k))
            #     hapmin[i] = score
            # end

            # keep all previous best pairs and equally good pairs
            # if score <= hapmin[i]
            #     push!(happairs[i], (j, k))
            #     hapmin[i] = score
            # end
        end
    end

    return nothing
end

"""
    fillmissing!(Xm, Xwork, H, haplopairs)

Fill in missing genotypes in `X` according to haplotypes. Non-missing genotypes
remain same.

# Input
* `Xm`: `p x n` genotype matrix with missing values. Each column is an individual.
* `Xwork`: `p x n` genotype matrix where missing values are filled with sum of 2 haplotypes.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `happair`: pair of haplotypes. `X[:, k] = H[:, happair[1][k]] + H[:, happair[2][k]]`.
"""
function fillmissing!(
    Xm::AbstractMatrix{Union{U, Missing}},
    Xwork::AbstractMatrix{T},
    H::AbstractMatrix{T},
    happairs::Vector{Vector{Tuple{Int, Int}}},
    ) where {T, U}

    p, n = size(Xm)
    best_discrepancy = typemax(eltype(Xwork))
    
    for j in 1:n, happair in happairs[j]
        discrepancy = zero(T)
        for i in 1:p
            if ismissing(Xm[i, j])
                tmp = H[i, happair[1]] + H[i, happair[2]]
                discrepancy += abs2(Xwork[i, j] - tmp)
                Xwork[i, j] = tmp
            end
        end
        if discrepancy < best_discrepancy
            best_discrepancy = discrepancy
        end
    end
    return best_discrepancy
end

# """
#     fillgeno!(X, H, happair)

# Fill in genotypes according to haplotypes. Both missing and non-missing
# genotypes may be changed.

# # Input
# * `X`: `p x n` genotype matrix. Each column is an individual.
# * `H`: `p x d` haplotype matrix. Each column is a haplotype.
# * `happair`: pair of haplotypes. `X[:, k] = H[:, happair[1][k]] + H[:, happair[2][k]]`.
# """
# function fillgeno!(
#     X::AbstractMatrix,
#     H::AbstractMatrix,
#     happair::Tuple{AbstractVector, AbstractVector}
#     )

#     @inbounds for j in 1:size(X, 2), i in 1:size(X, 1)
#         X[i, j] = H[i, happair[1][j]] + H[i, happair[2][j]]
#     end
#     return nothing

# end

"""
    initXfloat(X, Xfloat)

Initializes the matrix `Xfloat` where missing values of matrix `X` by `2 x` allele frequency
and nonmissing entries of `X` are converted to type `Float32` for subsequent BLAS routines. 

# Input
* `X` is a `p x n` genotype matrix. Each column is an individual.
* `Xfloat` is the `p x n` matrix of X where missing values are filled by 2x allele frequency. 
"""
function initXfloat!(
    X::AbstractMatrix;
    Xfloat::AbstractMatrix = zeros(Float32, size(X)),
    )
    
    T = Float32
    p, n = size(X)

    @inbounds for i in 1:p
        # allele frequency
        cnnz = zero(T)
        csum = zero(T)
        for j in 1:n
            if !ismissing(X[i, j])
                cnnz += one(T)
                csum += convert(T, X[i, j])
            end
        end
        # set missing values to 2freq, unless cnnz is 0
        imp = (cnnz == 0 ? zero(T) : csum / cnnz)
        for j in 1:n
            if ismissing(X[i, j])
                Xfloat[i, j] = imp
            else
                Xfloat[i, j] = convert(T, X[i, j])
            end
        end
    end

    any(isnan, Xfloat) && error("Xfloat contains NaN during initialization! Shouldn't happen!")
    any(isinf, Xfloat) && error("Xfloat contains Inf during initialization! Shouldn't happen!")
    any(ismissing, Xfloat) && error("Xfloat contains Missing during initialization! Shouldn't happen!")

    # impute using mode
    # for i in 1:p
    #     # set missing values to mode
    #     imp = mode(@view(X[i, :]))
    #     for j in 1:n
    #         if ismissing(X[i, j]) 
    #             Xfloat[i, j] = imp
    #         else
    #             Xfloat[i, j] = X[i, j]
    #         end
    #     end
    # end

    # initialize using 0
    # for i in 1:p, j in 1:n
    #     Xfloat[i, j] = ifelse(ismissing(X[i, j]), zero(T), X[i, j])
    # end

    return nothing
end

"""
    haploimpute!(X, H, M, N, happair, hapscore, maxiters=1, tolfun=1e-3)

In a window, performs haplotying of genotype matrix `X` from the pool of 
haplotypes `H`.

# Input
* `X`: `p x n` matrix with missing values. Each column is genotypes of an individual.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `M`: overwritten by `M[i, j] = 2dot(H[:, i], H[:, j]) + sumabs2(H[:, i]) +
    sumabs2(H[:, j])`.
* `N`: overwritten by `n x d` matrix `2X'H`.
* `happair`: vector of optimal haplotype pair. `X[:, k] ≈ H[:, happair[k, 1]] + H[:, happair[k, 2]]`.
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
* `Xfloat`: copy of `X` where missing values are filled with mean. This engages in linear algebra for computing `N`
* `maxiters`: number of MM iterations. Default is 1.
* `tolfun`: convergence tolerance of MM iterations. Default is 1e-3.
"""
function haploimpute!(
    X::AbstractMatrix,
    H::AbstractMatrix{T},
    M::AbstractMatrix{T},
    N::AbstractMatrix{T},
    happairs::Vector{Vector{Tuple{Int, Int}}},
    hapscore::AbstractVector;
    Xfloat::AbstractMatrix = zeros(T, size(X)),
    maxiters::Int  = 1,
    tolfun::Number = 1e-3,
    ) where T

    obj = typemax(eltype(hapscore))
    size(X) == size(Xfloat) || error("Dimension mismatch: X and Xfloat have sizes $(size(X)) and $(size(Xfloat))")
    initXfloat!(X, Xfloat=Xfloat) #Xfloat is the matrix that engages in BLAS routines

    # mm iteration
    for iter in 1:maxiters
        # compute top haplotype pairs for each genotype vector
        haplopair!(Xfloat, H, M, N, happairs, hapscore)
        # screen for best haplotype pair based on observed entries
        choose_happair!(X, H, happairs, hapscore)
        # impute missing entries according to current haplotypes
        discrepancy = fillmissing!(X, Xfloat, H, happairs)
        # convergence criterion
        objold = obj
        obj = sum(hapscore) - discrepancy
        # println("iter = $iter, discrepancy = $discrepancy, obj = $obj")
        if abs(obj - objold) < tolfun * (objold + 1)
            break
        end
    end
    return nothing
end

"""
    choose_happair!(X, H, happairs, hapscore)

Calculates error ||x - hi - hj||^2 only on the observed entries and save 
observed error in `hapscore`. `happairs` will keep only the best haplotype 
pairs based on the error of observed entries. All happairs
that attain the best observed error will be kept.
"""
function choose_happair!(
    X::AbstractMatrix{Union{Missing, T}},
    H::AbstractMatrix,
    happairs::Vector{Vector{Tuple{Int, Int}}},
    hapscore::AbstractVector;
    ) where T <: Real

    p = size(X, 1)
    n = size(X, 2)
    d = size(H, 2)
    p == size(H, 1) || error("Dimension mismatch: size(X, 1) = $p but size(H, 1) = $(size(H, 1))")
    U = promote_type(eltype(X), eltype(H))

    # loop over each person's genotype
    best_happair = Tuple{Int, Int}[]
    for j in 1:n
        best_error = typemax(eltype(H))
        empty!(best_happair)
        for happair in happairs[j]
            # compute errors for each pair based on observed entries
            h1, h2 = happair[1], happair[2]
            err = zero(U)
            @inbounds @simd for i in 1:p
                if X[i, j] !== missing 
                    err += (X[i, j] - H[i, h1] - H[i, h2])^2
                end
            end
            if err :: U < best_error
                best_error = err
                empty!(best_happair)
                push!(best_happair, happair)
            elseif err :: U == best_error
                push!(best_happair, happair)
            end
        end

        # keep only best haplotype pair in happairs
        if length(happairs[j]) > 1
            empty!(happairs[j])
            for pair in best_happair
                push!(happairs[j], pair)
            end
        end
        hapscore[j] = convert(eltype(hapscore), best_error)
    end

    return nothing
end

"""
    chunk_size(people, haplotypes)

Figures out how many SNPs can be loaded into memory (capped at 2/3 total RAM) 
at once, given the data size. Assumes genotype data are Float32 (4 byte per entry) 
and haplotype panels are BitArrays (1 bit per entry).
"""
function chunk_size(people::Int, haplotypes::Int)
    system_memory_gb = Sys.total_memory() / 2^30
    system_memory_bits = 8000000000 * system_memory_gb
    usable_bits = round(Int, system_memory_bits * 2 / 3) # use 2/3 of memory for genotype and haplotype matrix per chunk
    max_chunk_size = round(Int, usable_bits / (haplotypes + 32people))
    return max_chunk_size

    # return 1000 # for testing in compare1
end
