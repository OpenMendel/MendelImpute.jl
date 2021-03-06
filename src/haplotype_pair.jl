###### This file is part of the MendelImpute.jl package.
###### It contains relevant code to solve the least squares criteria 
###### `x ≈ hᵢ + hⱼ`. 

"""
    compute_optimal_haplotypes!(haplotype1, haplotype2, compressed_Hunique, ...)

Computes the best haplotype pair `(hᵢ, hⱼ)` for each genotype vector `x` and 
stores result in `haplotype1` and `haplotype2`.

# Arguments
- `haplotype1`: Person `i` strand1 haplotype in window `w` is `haplotype1[i][w]`
- `haplotype2`: Person `i` strand2 haplotype in window `w` is `haplotype2[i][w]`
- `compressed_Hunique`: A `CompressedHaplotypes` object
- `X`: the full genotype matrix possibly with missings. Each column is an 
    individual.
- `X_pos`: Position of each SNP in `X`
- `typed_snps_error`: `typed_snps_error[i]` is the % of SNPs correctly imputed.
    `1` is best, `0` is worse.
- `stepscreen`: Boolean indicating whether to use a stepwise heuristic to screen
    for top haplotypes instead of performing global search
- `tf`: This option solves the least squares objective on only
    `tf` unique haplotypes.
- `scale_allelefreq` Boolean indicating whether to give rare SNPs more weight
    scaled by `wᵢ = 1 / √2p(1-p)` where max weight is 2. 
- `max_haplotypes` Maximum number of haplotypes for using to global search. 
    This number should be specified along with `stepscreen` or `tf`.
- `rescreen` Boolean for a more rigorous global search. 
- `timers`

# Timers:
- `t1` = screening for top haplotypes (0 unless `tf` option is on)
- `t2` = BLAS3 mul! to get M and N
- `t3` = haplopair search
- `t4` = rescreen time (0 unless `rescreen = true`)
- `t5` = initializing missing
- `t6` = allocating and viewing
- `t7` = index conversion
"""
function compute_optimal_haplotypes!(
    haplotype1::AbstractVector,
    haplotype2::AbstractVector,
    compressed_Hunique::CompressedHaplotypes,
    X::AbstractMatrix,
    X_pos::AbstractVector,
    stepscreen::Union{Nothing, Int},
    tf::Union{Nothing, Int}, # thinning factor
    scale_allelefreq::Bool,
    max_haplotypes::Int,
    rescreen::Bool
    )
    # constants
    people = size(X, 2)
    ref_snps = length(compressed_Hunique.pos)
    max_width, max_d = MendelImpute.max_dim(compressed_Hunique)
    windows = length(haplotype1[1])
    threads = Threads.nthreads()
    inv_sqrt_allele_var = nothing
    winranges = compressed_Hunique.X_window_range
    overlap = compressed_Hunique.overlap

    # allocate working arrays
    timers = [zeros(7*8) for _ in 1:threads] # 8 for spacing
    timers[1][48] += @elapsed begin # time for allocating
        happair1 = [ones(Int32, people)               for _ in 1:threads]
        happair2 = [ones(Int32, people)               for _ in 1:threads]
        hapscore = [zeros(Float32, people)            for _ in 1:threads]
        Xwork    = [zeros(Float32, max_width, people) for _ in 1:threads]
        Hwork    = [zeros(Float32, max_width, max_d)  for _ in 1:threads]
        if !isnothing(tf)
            maxindx = [zeros(Int32, tf)         for _ in 1:threads]
            maxgrad = [zeros(Float32, tf)       for _ in 1:threads]
            Hk = [zeros(Float32, max_width, tf) for _ in 1:threads]
            Xi = [zeros(Float32, max_width)     for _ in 1:threads]
            M  = [zeros(Float32, tf, tf)        for _ in 1:threads]
            N  = [zeros(Float32, tf)            for _ in 1:threads]
            R  = [zeros(Float32, max_d, people) for _ in 1:threads]
        elseif !isnothing(stepscreen)
            maxindx = [zeros(Int32,   stepscreen) for _ in 1:threads]
            maxgrad = [zeros(Float32, stepscreen) for _ in 1:threads]
            M  = [zeros(Float32, max_d, max_d)    for _ in 1:threads]
            Nt = [zeros(Float32, max_d, people)   for _ in 1:threads]
        else
            M = [zeros(Float32, max_d,  max_d) for _ in 1:threads]
            N = [zeros(Float32, people, max_d) for _ in 1:threads]
        end
    end
    pmeter = Progress(windows, 5, "Computing optimal haplotypes...")

    # for w in 1:windows
    Threads.@threads for w in 1:windows
        id = Threads.threadid()
        t6 = @elapsed begin
            Hw_aligned = compressed_Hunique.CW_typed[w].uniqueH
            Xrange = extend_to_overlap_range(compressed_Hunique, w, overlap)
            Xw_aligned = view(X, Xrange, :)
            d = size(Hw_aligned, 2)
        end

        # weight snp by inverse allele variance if requested
        t6 += @elapsed if scale_allelefreq
            Hw_range = compressed_Hunique.Hstart[w]:(w ==
                windows ? ref_snps : compressed_Hunique.Hstart[w + 1] - 1)
            Hw_snp_pos = indexin(X_pos[winranges[w]],
                compressed_Hunique.pos[Hw_range])
            inv_sqrt_allele_var = compressed_Hunique.altfreq[Hw_snp_pos]
            map!(x -> x < 0.15 ? 1.98 : 1 / sqrt(2*x*(1-x)),
                inv_sqrt_allele_var, inv_sqrt_allele_var) # set min pᵢ = 0.15
        end

        # compute top haplotype pairs for each sample in current window
        if !isnothing(stepscreen) && d > max_haplotypes
            # find hᵢ via stepwise regression, then find hⱼ via global search
            t1, t2, t3, t4, t5, t6 = haplopair_stepscreen!(Xw_aligned, 
                Hw_aligned, r=stepscreen, 
                inv_sqrt_allele_var=inv_sqrt_allele_var, happair1=happair1[id], 
                happair2=happair2[id], hapscore=hapscore[id], 
                maxindx=maxindx[id], maxgrad=maxgrad[id], Xwork=Xwork[id],
                Hwork=Hwork[id], Nt=Nt[id], M=M[id])
        elseif !isnothing(tf) && d > max_haplotypes
            # haplotype thinning: search all (hᵢ, hⱼ) pairs where hᵢ ≈ x ≈ hⱼ
            t1, t2, t3, t4, t5, t6 = haplopair_thin_BLAS2!(Xw_aligned,
                Hw_aligned, allele_freq=inv_sqrt_allele_var, keep=tf,
                happair1=happair1[id], happair2=happair2[id],
                hapscore=hapscore[id], maxindx=maxindx[id], maxgrad=maxgrad[id],
                Xi=Xi[id], N=N[id], Hk=Hk[id], M=M[id], Xwork=Xwork[id],
                Hwork=Hwork[id], R=R[id])
        elseif rescreen
            # global search + searching ||x - hᵢ - hⱼ|| on observed entries
            t1, t2, t3, t4, t5, t6 = haplopair_rescreen!(Xw_aligned, 
                Hw_aligned, happair1=happair1[id], happair2=happair2[id],
                hapscore=hapscore[id], Xwork=Xwork[id])
        else
            # global search
            t1, t2, t3, t4, t5, t6 = haplopair!(Xw_aligned, Hw_aligned,
                inv_sqrt_allele_var=inv_sqrt_allele_var, happair1=happair1[id],
                happair2=happair2[id], hapscore=hapscore[id], Xwork=Xwork[id],
                Hwork=Hwork[id], M=M[id], N = N[id])
        end

        # save result 
        t7 = @elapsed save_haplotypes!(haplotype1, haplotype2, happair1[id], 
            happair2[id], compressed_Hunique, w)

        # record timings and haplotypes (× 8 to avoid false sharing)
        timers[id][8]  += t1
        timers[id][16] += t2
        timers[id][24] += t3
        timers[id][32] += t4
        timers[id][40] += t5
        timers[id][48] += t6
        timers[id][56] += t7

        # update progress
        next!(pmeter)
    end

    return sum(timers) ./ threads
end

"""
    save_haplotypes!(haplotype1, haplotype2, happair1, happair2, ...)

Helper function to convert `happair`s (which index off unique haplotypes) to 
indices of full haplotype pool, and store them in `haplotype`s.

# Arguments
- `haplotype1`: Person `i` strand1 haplotype in window `w` is `haplotype1[i][w]`
- `haplotype2`: Person `i` strand2 haplotype in window `w` is `haplotype2[i][w]`
- `happair1`: Optimal haplotype pair in strand1 of current window, indexes off
    of unique haplotypes.
- `happair2`: Optimal haplotype pair in strand2 of current window, indexes off
    of unique haplotypes.
- `compressed_Hunique`: A `CompressedHaplotypes` object
- `window` current window.
"""
function save_haplotypes!(
    haplotype1::Vector{Vector{Int32}},
    haplotype2::Vector{Vector{Int32}},
    happair1::Vector{Int32},
    happair2::Vector{Int32},
    compressed_Hunique::CompressedHaplotypes,
    window::Int,
    )

    people = length(haplotype1)
    @inbounds for i in 1:people
        haplotype1[i][window] = unique_idx_to_complete_idx(
            happair1[i], window, compressed_Hunique)
        haplotype2[i][window] = unique_idx_to_complete_idx(
            happair2[i], window, compressed_Hunique)
    end
    return nothing
end

# """
# Records optimal-redundant haplotypes for each window.

# Warning: This function is called in a multithreaded loop. If you modify this
# function you must check whether imputation accuracy is affected (when run with
# >1 threads).

# # Arguments:
# - `window_idx`: window in current chunk
# - `window_overall`: window index in terms of every windows
# """
# function compute_redundant_haplotypes!(
#     redundant_haplotypes::Vector{OptimalHaplotypeSet},
#     Hunique::CompressedHaplotypes,
#     happair1::AbstractVector,
#     happair2::AbstractVector,
#     window_idx::Int,
#     window_overall::Int,
#     storage1 = falses(nhaplotypes(Hunique)),
#     storage2 = falses(nhaplotypes(Hunique))
#     )

#     people = length(redundant_haplotypes)

#     @inbounds for k in 1:people
#         # convert happairs from unique idx to complete idx
#         Hi_idx = unique_idx_to_complete_idx(happair1[k], window_overall,
#             Hunique)
#         Hj_idx = unique_idx_to_complete_idx(happair2[k], window_overall,
#             Hunique)

#         # strand1
#         storage1 .= false
#         if haskey(Hunique.CW_typed[window_overall].hapmap, Hi_idx)
#             h1_set = Hunique.CW_typed[window_overall].hapmap[Hi_idx]
#             for i in h1_set
#                 storage1[i] = true
#             end
#         else
#             storage1[Hi_idx] = true # Hi_idx is singleton (i.e. unique)
#         end

#         # strand2
#         storage2 .= false
#         if haskey(Hunique.CW_typed[window_overall].hapmap, Hj_idx)
#             h2_set = Hunique.CW_typed[window_overall].hapmap[Hj_idx]
#             for i in h2_set
#                 storage2[i] = true
#             end
#         else
#             storage2[Hj_idx] = true # Hj_idx is singleton (i.e. unique)
#         end

#         # redundant_haplotypes[k].strand1[window_idx] = copy(storage1)
#         # redundant_haplotypes[k].strand2[window_idx] = copy(storage2)
#         if isassigned(redundant_haplotypes[k].strand1, window_idx)
#             redundant_haplotypes[k].strand1[window_idx] .= storage1
#             redundant_haplotypes[k].strand2[window_idx] .= storage2
#         else
#             redundant_haplotypes[k].strand1[window_idx] = copy(storage1)
#             redundant_haplotypes[k].strand2[window_idx] = copy(storage2)
#         end
#     end

#     return nothing
# end

"""
    screen_flanking_windows!(haplotype1, haplotype2, compressed_Hunique, X)

For each window's haplotype pair, tests whether adjacent window's haplotype
pairs produce better error on the observed entries. 

# Arguments
- `haplotype1`: Person `i` strand1 haplotype in window `w` is `haplotype1[i][w]`
- `haplotype2`: Person `i` strand2 haplotype in window `w` is `haplotype2[i][w]`
- `compressed_Hunique`: A `CompressedHaplotypes` object
- `X`: the full genotype matrix possibly with missings. Each column is an
    individual.
"""
function screen_flanking_windows!(
    haplotype1::AbstractVector,
    haplotype2::AbstractVector,
    compressed_Hunique::CompressedHaplotypes,
    X::AbstractMatrix
    )

    people = length(haplotype1)
    haplotypes = nhaplotypes(compressed_Hunique)
    width = compressed_Hunique.width
    windows = length(haplotype1[1])
    winranges = compressed_Hunique.X_window_range

    for w in 1:windows
        Hw_aligned = compressed_Hunique.CW_typed[w].uniqueH
        Xw_aligned = view(X, winranges[w], :)

        for i in 1:people
            # calculate observed error for current pair
            h1_curr_complete = haplotype1[i][w] # complete index
            h2_curr_complete = haplotype2[i][w] # complete index
            h1_curr = complete_idx_to_unique_typed_idx(h1_curr_complete, w, 
                compressed_Hunique) # unique index
            h2_curr = complete_idx_to_unique_typed_idx(h2_curr_complete, w, 
                compressed_Hunique) # unique index
            curr_err = observed_error(Xw_aligned, i, Hw_aligned, h1_curr, 
                h2_curr)

            # consider previous pair
            if w != 1
                h1_prev_complete = haplotype1[i][w - 1]
                h2_prev_complete = haplotype2[i][w - 1]
                h1_prev = complete_idx_to_unique_typed_idx(h1_prev_complete, 
                    w, compressed_Hunique)
                h2_prev = complete_idx_to_unique_typed_idx(h2_prev_complete, 
                    w, compressed_Hunique)
                prev_err = observed_error(Xw_aligned, i, Hw_aligned, h1_prev, 
                    h2_prev)
                if prev_err < curr_err
                    h1_curr, h2_curr, curr_err = h1_prev, h2_prev, prev_err
                end
            end

            # consider next pair
            if w != windows
                h1_next_complete = haplotype1[i][w + 1]
                h2_next_complete = haplotype2[i][w + 1]
                h1_next = complete_idx_to_unique_typed_idx(h1_next_complete, 
                    w, compressed_Hunique) # unique index
                h2_next = complete_idx_to_unique_typed_idx(h2_next_complete, 
                    w, compressed_Hunique) # unique index
                next_err = observed_error(Xw_aligned, i, Hw_aligned, h1_next, 
                    h2_next)
                if next_err < curr_err
                    h1_curr, h2_curr, curr_err = h1_next, h2_next, next_err
                end
            end

            # convert from unique idx back to complete idx
            H1_idx = unique_idx_to_complete_idx(h1_curr, w, compressed_Hunique)
            H2_idx = unique_idx_to_complete_idx(h2_curr, w, compressed_Hunique)
            haplotype1[i][w] = H1_idx
            haplotype2[i][w] = H2_idx
        end
    end

    return nothing
end

# uses dynamic programming. Only the first 1000 haplotype pairs will be saved.
# function compute_redundant_haplotypes!(
#     redundant_haplotypes::Vector{Vector{Vector{T}}},
#     Hunique::CompressedHaplotypes,
#     happair1::AbstractVector,
#     happair2::AbstractVector,
#     window_idx::Int,
#     window_overall::Int,
#     storage1 = falses(nhaplotypes(Hunique)),
#     storage2 = falses(nhaplotypes(Hunique))
#     ) where T <: Tuple{Int32, Int32}

#     people = length(redundant_haplotypes)

#     @inbounds for k in 1:people
#         # convert happairs from unique idx to complete idx
#         Hi_idx = unique_idx_to_complete_idx(happair1[k], window_overall, Hunique)
#         Hj_idx = unique_idx_to_complete_idx(happair2[k], window_overall, Hunique)

#         # find haplotypes that match Hi_idx and Hj_idx on typed snps
#         h1_set = get(Hunique.CW_typed[window_overall].hapmap, Hi_idx, Hi_idx)
#         h2_set = get(Hunique.CW_typed[window_overall].hapmap, Hj_idx, Hj_idx)

#         # save first 1000 haplotype pairs
#         for h1 in h1_set, h2 in h2_set
#             if length(redundant_haplotypes[k][window_idx]) < 1000
#                 push!(redundant_haplotypes[k][window_idx], (h1, h2))
#             else
#                 break
#             end
#         end
#     end

#     return nothing
# end

"""
    haplopair(X, H)

Calculate the best pair of haplotypes in `H` for each individual in `X`. 
Missing data is initialized as 2x alternate allele freq.

# Input
* `X`: `p × n` genotype matrix possibly with missings. Each column is an individual.
* `H`: `p × d` haplotype matrix. Each column is a haplotype.

# Optional arguments
* `happair1`: optimal haplotype pairs for strand 1.
* `happair2`: optimal haplotype pairs for strand 2.
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
* `inv_sqrt_allele_var`: If not nothing, SNP `i` will be scaled according to 
    `inv_sqrt_allele_var[i]`.
* `Xwork`: A `Float32` version of `X`, missing entries are filled with `2maf`
* `Hwork`: A `Float32` version of `H`.
* `M`: `d x d` matrix with entries `M[i, j] = 2dot(Hwork[:, i], Hwork[:, j]) +
    sumabs2(Hwork[:, i]) + sumabs2(Hwork[:, j])`. Only the upper triangular 
    part of `M` is used.
* `N`: `n x d` matrix `2Xwork'Hwork`.
"""
function haplopair!(
    X::AbstractMatrix, # p × n
    H::AbstractMatrix; # p × d
    # preallocated vectors
    happair1::AbstractVector = ones(Int32, size(X, 2)), # length n
    happair2::AbstractVector = ones(Int32, size(X, 2)), # length n
    hapscore::AbstractVector = Vector{Float32}(undef, size(X, 2)), # length n
    inv_sqrt_allele_var::Union{Nothing, AbstractVector} = nothing, # length p
    # preallocated matrices
    Xwork::AbstractMatrix{Float32}=Matrix{Float32}(undef,size(X,1),size(X,2)), # p × n
    Hwork::AbstractMatrix{Float32}=convert(Matrix{Float32}, H),                # p × d
    M::AbstractMatrix{Float32}=Matrix{Float32}(undef, size(H, 2), size(H, 2)), # d × d
    N::AbstractMatrix{Float32}=Matrix{Float32}(undef, size(X, 2), size(H, 2)), # n × d
    )
    p,  n = size(X)
    pp, d = size(H)
    p == pp || error("size of H doesn't match X! Expected $p got $pp")

    # create views for matrices
    t6 = @elapsed begin
        Mwork_view = view(M, 1:d, 1:d)
        Nwork_view = view(N, :, 1:d)
        Hwork_view = view(Hwork, 1:p, 1:d)
        Xwork_view = view(Xwork, 1:p, :)
        copyto!(Hwork_view, H) # update Hwork_view, others are updated later
    end

    # initializes missing
    t5 = @elapsed initXfloat!(Xwork_view, X)

    t2, t3 = haplopair!(Xwork_view, Hwork_view, Mwork_view, Nwork_view, 
        happair1, happair2, hapscore, inv_sqrt_allele_var)
    t1 = t4 = 0.0 # no time spent on haplotype thinning or rescreening

    return t1, t2, t3, t4, t5, t6
end

"""
    haplopair!(X, H, M, N, happair, hapscore)

Calculate the best pair of haplotypes in `H` for each individual in `X`.
Overwite `M` by `M[i, j] = 2dot(H[:, i], H[:, j]) + sumabs2(H[:, i]) + 
sumabs2(H[:, j])`, `N` by `2X'H`, `happair` by optimal haplotype pair, and 
`hapscore` by objective value from the optimal haplotype pair.

# Input
* `X`: `p x n` genotype matrix. Each column is an individual.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `M`: overwritten by `M[i, j] = 2dot(H[:, i], H[:, j]) + sumabs2(H[:, i]) +
    sumabs2(H[:, j])`.
* `N`: overwritten by `n x d` matrix `2X'H`.
* `happair`: optimal haplotype pair. `X[:, k] ≈ H[:, happair[k, 1]] + 
    H[:, happair[k, 2]]`.
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
"""
function haplopair!(
    X::AbstractMatrix{Float32},
    H::AbstractMatrix{Float32},
    M::AbstractMatrix{Float32},
    N::AbstractMatrix{Float32},
    happair1::AbstractVector{Int32},
    happair2::AbstractVector{Int32},
    hapscore::AbstractVector{Float32},
    inv_sqrt_allele_var::Union{Nothing, AbstractVector}
    )

    p, n, d = size(X, 1), size(X, 2), size(H, 2)

    # assemble M (upper triangular only)
    t2 = @elapsed begin
        if !isnothing(inv_sqrt_allele_var)
            H .*= inv_sqrt_allele_var # wᵢ = 1/√2p(1-p)
        end
        mul!(M, Transpose(H), H)
        for j in 1:d, i in 1:(j - 1) # off-diagonal
            @inbounds M[i, j] = 2M[i, j] + M[i, i] + M[j, j]
        end
        for j in 1:d # diagonal
            @inbounds M[j, j] *= 4
        end

        # assemble N
        if !isnothing(inv_sqrt_allele_var)
            H .*= inv_sqrt_allele_var # wᵢ = 1/2p(1-p)
        end
        mul!(N, Transpose(X), H)
        @simd for I in eachindex(N)
            N[I] *= 2
        end
    end

    # computational routine
    t3 = @elapsed haplopair!(happair1, happair2, hapscore, M, N)

    # supplement the constant terms in objective
    t3 += @elapsed begin
        @inbounds for j in 1:n
            @simd for i in 1:p
                hapscore[j] += abs2(X[i, j])
            end
        end
    end

    return t2, t3
end

"""
    haplopair!(happair, hapscore, M, N)

Calculate the best pair of haplotypes pairs in the filtered haplotype panel
for each individual in `X` using sufficient statistics `M` and `N`.

# Note
The best haplotype pairs are column indices of the filtered haplotype panels.

# Input
* `happair1`: optimal haplotype in strand 1 for each individual.
* `happair2`: optimal haplotype in strand 2 for each individual.
* `hapmin`: minimum offered by the optimal haplotype pair.
* `M`: `d x d` matrix with entries `M[i, j] = 2dot(H[:, i], H[:, j]) +
    sumabs2(H[:, i]) + sumabs2(H[:, j])`, where `H` is the haplotype matrix
    with haplotypes in columns. Only the upper triangular part of `M` is used.
* `N`: `n x d` matrix `2X'H`, where `X` is the genotype matrix with individuals
    in columns.
"""
function haplopair!(
    happair1::AbstractVector{Int32},
    happair2::AbstractVector{Int32},
    hapmin::AbstractVector{Float32},
    M::AbstractMatrix{Float32},
    N::AbstractMatrix{Float32},
    )

    n, d = size(N)
    chunks = div(n, 4)
    fill!(hapmin, Inf32)
    @inbounds for k in 1:d, j in 1:k
        Mjk = M[j, k]
        i = 1
        # loop over individuals
        for chunk in 1:chunks # unrolled loop, equivalent to for i in 1:n
            i2, i3, i4 = i + 1, i + 2, i + 3
            score = Mjk - N[i, j] - N[i, k]
            if score < hapmin[i]
                hapmin[i], happair1[i], happair2[i] = score, j, k
            end
            score = Mjk - N[i2, j] - N[i2, k]
            if score < hapmin[i2]
                hapmin[i2], happair1[i2], happair2[i2] = score, j, k
            end
            score = Mjk - N[i3, j] - N[i3, k]
            if score < hapmin[i3]
                hapmin[i3], happair1[i3], happair2[i3] = score, j, k
            end
            score = Mjk - N[i4, j] - N[i4, k]
            if score < hapmin[i4]
                hapmin[i4], happair1[i4], happair2[i4] = score, j, k
            end
            i += 4
        end
        # handle remaining terms
        while i ≤ n
            score = Mjk - N[i, j] - N[i, k]
            if score < hapmin[i]
                hapmin[i], happair1[i], happair2[i] = score, j, k
            end
            i += 1
        end
    end
end

"""
    fillmissing!(Xm, Xwork, H, haplopairs)

Fill in missing genotypes in `X` according to haplotypes. Non-missing genotypes
remain same.

# Input
* `Xm`: `p x n` genotype matrix with missing values. Each column is an
    individual.
* `Xwork`: `p x n` genotype matrix where missing values are filled with sum of 2
    haplotypes.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `happair`: pair of haplotypes. `X[:, k] = H[:, happair[1][k]] + 
    H[:, happair[2][k]]`.
"""
function fillmissing!(
    Xm::AbstractMatrix{Union{U, Missing}},
    Xwork::AbstractMatrix{T},
    H::AbstractMatrix{T},
    happairs::Vector{Vector{Tuple{Int32, Int32}}},
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

"""
    initXfloat!(Xfloat, X)

Initializes the matrix `Xfloat` where missing values of matrix `X` by `2 x`
allele frequency and nonmissing entries of `X` are converted to type `Float32`
for subsequent BLAS routines.

# Input
* `X` is a `p x n` genotype matrix. Each column is an individual.
* `Xfloat` is the `p x n` matrix of X where missing values are filled by 2x
    allele frequency.
"""
function initXfloat!(
    Xfloat::AbstractMatrix,
    X::AbstractMatrix
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

    any(isnan, Xfloat) && error("Xfloat has NaN during initialization!")
    any(isinf, Xfloat) && error("Xfloat has Inf during initialization!")
    any(ismissing, Xfloat) && error("Xfloat has Missing during initialization!")

    return nothing
end
