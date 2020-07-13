"""
    haplochunk!(redundant_haplotypes, compressed_Hunique, X, ...)

Wrapper function that computes the best haplotype pair `(hᵢ, hⱼ)` for each genotype vector
in a given chunk.

# There are 5 timers (some may be 0):
t1 = screening for top haplotypes
t2 = BLAS3 mul! to get M and N
t3 = haplopair search
t4 = rescreen time
t5 = finding redundant happairs
"""
function haplochunk!(
    redundant_haplotypes::AbstractVector,
    compressed_Hunique::CompressedHaplotypes,
    X::AbstractMatrix,
    X_pos::AbstractVector,
    dynamic_programming::Bool,
    lasso::Union{Nothing, Int},
    thinning_factor::Union{Nothing, Int},
    thinning_scale_allelefreq::Bool,
    max_haplotypes::Int,
    rescreen::Bool,
    winrange::UnitRange,
    total_window::Int,
    pmeter::Progress,
    timers::AbstractVector
    )
    people = size(X, 2)
    ref_snps = length(compressed_Hunique.pos)
    width = compressed_Hunique.width
    windows = length(winrange)
    threads = Threads.nthreads()
    avghaps = nhaplotypes(compressed_Hunique)

    # working arrys 
    happair1 = [ones(Int32, people)           for _ in 1:threads]
    happair2 = [ones(Int32, people)           for _ in 1:threads]
    hapscore = [zeros(Float32, size(X, 2))    for _ in 1:threads]
    Xwork    = [zeros(Float32, width, people) for _ in 1:threads]
    # N        = [ElasticArray{Float32}(undef, people, avghaps) for _ in 1:threads]
    # Hwork    = [ElasticArray{Float32}(undef, width, avghaps)  for _ in 1:threads]
    if !isnothing(thinning_factor)
        maxindx = [zeros(Int, thinning_factor)                 for _ in 1:threads]
        maxgrad = [zeros(Float32, thinning_factor)             for _ in 1:threads]
        Hk = [zeros(Float32, width, thinning_factor)           for _ in 1:threads]
        Xi = [zeros(Float32, width)                            for _ in 1:threads]
        M  = [zeros(Float32, thinning_factor, thinning_factor) for _ in 1:threads]
        N  = [zeros(Float32, thinning_factor)                  for _ in 1:threads]
    end

    ThreadPools.@qthreads for absolute_w in winrange
        Hw_aligned = compressed_Hunique.CW_typed[absolute_w].uniqueH
        Xw_idx_start = (absolute_w - 1) * width + 1
        Xw_idx_end = (absolute_w == total_window ? length(X_pos) : absolute_w * width)
        Xw_aligned = view(X, Xw_idx_start:Xw_idx_end, :)
        id = Threads.threadid()

        # computational routine
        if !isnothing(lasso)
            if r > size(Hw_aligned, 2) || size(Hw_aligned, 2) <= max_haplotypes
                # global search
                t1, t2, t3, t4 = haplopair!(Xw_aligned, Hw_aligned, happair1=happair1[id], 
                    happair2=happair2[id], hapscore=hapscore[id])
            else
                happairs, hapscore, t1, t2, t3, t4 = haplopair_lasso(Xw_aligned, Hw_aligned, 
                    r = lasso)
            end
        elseif !isnothing(thinning_factor)
            # weight each snp by frequecy if requested
            if thinning_scale_allelefreq
                Hw_range = compressed_Hunique.start[absolute_w]:(absolute_w == total_window ? 
                    ref_snps : compressed_Hunique.start[absolute_w + 1] - 1)
                Hw_snp_pos = indexin(X_pos[Xw_idx_start:Xw_idx_end], compressed_Hunique.pos[Hw_range])
                altfreq = compressed_Hunique.altfreq[Hw_snp_pos]
            else
                altfreq = nothing
            end
            # run haplotype thinning (i.e. search all (hi, hj) pairs where hi, hj ≈ x)
            if thinning_factor < size(Hw_aligned, 2) || size(Hw_aligned, 2) > max_haplotypes
                t1, t2, t3, t4 = haplopair_thin_BLAS2!(Xw_aligned, Hw_aligned, 
                    alt_allele_freq=altfreq, keep=thinning_factor, happair1=happair1[id], 
                    happair2=happair2[id], hapscore=hapscore[id], maxindx=maxindx[id], 
                    maxgrad=maxgrad[id], Xi=Xi[id], N=N[id], Hk=Hk[id], M=M[id], 
                    Xwork=Xwork[id])
            else
                t1, t2, t3, t4 = haplopair!(Xw_aligned, Hw_aligned, happair1=happair1[id], 
                    happair2=happair2[id], hapscore=hapscore[id])
            end
        elseif rescreen
            # global search to find many (hi, hj) pairs, then reminimize ||x - hi - hj|| on observed entries
            happairs, hapscore, t1, t2, t3, t4 = haplopair_screen(Xw_aligned, Hw_aligned)
        else
            # global search
            t1, t2, t3, t4 = haplopair!(Xw_aligned, Hw_aligned, happair1=happair1[id], 
                happair2=happair2[id], hapscore=hapscore[id], Xwork=Xwork[id])
        end

        # convert happairs (which index off unique haplotypes) to indices of full haplotype pool, and find all matching happairs
        t5 = @elapsed begin
            w = something(findfirst(x -> x == absolute_w, winrange)) # window index of current chunk
            compute_redundant_haplotypes!(redundant_haplotypes, compressed_Hunique, happair1[id], 
                happair2[id], w, absolute_w, dp = dynamic_programming)
        end
         
        # record timings and haplotypes
        id = Threads.threadid()
        timers[id][1] += t1
        timers[id][2] += t2
        timers[id][3] += t3
        timers[id][4] += t4
        timers[id][5] += t5

        # update progress
        next!(pmeter)
    end
end

"""
Records optimal-redundant haplotypes for each window. Currently, only the first 1000
haplotype pairs will be saved to reduce search space for dynamic programming. 

Warning: This function is called in a multithreaded loop. If you modify this function
you must check whether imputation accuracy is affected (when run with >1 threads).

# Arguments:
- `window_idx`: window in current chunk
- `window_overall`: window index in terms of every windows
"""
function compute_redundant_haplotypes!(
    redundant_haplotypes::Union{Vector{Vector{Vector{T}}}, Vector{OptimalHaplotypeSet}}, 
    Hunique::CompressedHaplotypes, 
    happair1::AbstractVector,
    happair2::AbstractVector,
    window_idx::Int,
    window_overall::Int;
    dp::Bool = false, # dynamic programming
    ) where T <: Tuple{Int32, Int32}
    
    people = length(redundant_haplotypes)
    haplotypes = nhaplotypes(Hunique)

    if dp
        @inbounds for k in 1:people
            Hi_idx = unique_idx_to_complete_idx(happair1[k], window_overall, Hunique)
            Hj_idx = unique_idx_to_complete_idx(happair2[k], window_overall, Hunique)

            # find haplotypes that match Hi_idx and Hj_idx on typed snps
            h1_set = get(Hunique.CW_typed[window_overall].hapmap, Hi_idx, Hi_idx)
            h2_set = get(Hunique.CW_typed[window_overall].hapmap, Hj_idx, Hj_idx)

            # save first 1000 haplotype pairs
            for h1 in h1_set, h2 in h2_set
                if length(redundant_haplotypes[k][window_idx]) < 1000 
                    push!(redundant_haplotypes[k][window_idx], (h1, h2))
                else
                    break
                end
            end
        end
    else
        redunhaps_bitvec1 = falses(haplotypes)
        redunhaps_bitvec2 = falses(haplotypes)
        @inbounds for k in 1:people
            Hi_idx = unique_idx_to_complete_idx(happair1[k], window_overall, Hunique)
            Hj_idx = unique_idx_to_complete_idx(happair2[k], window_overall, Hunique)

            # find haplotypes that match Hi_idx and Hj_idx on typed snps
            h1_set = get(Hunique.CW_typed[window_overall].hapmap, Hi_idx, Hi_idx)
            h2_set = get(Hunique.CW_typed[window_overall].hapmap, Hj_idx, Hj_idx)

            # record matching haplotypes into bitvector
            redunhaps_bitvec1 .= false
            redunhaps_bitvec2 .= false
            for i in h1_set
                redunhaps_bitvec1[i] = true
            end
            for i in h2_set
                redunhaps_bitvec2[i] = true
            end
            redundant_haplotypes[k].strand1[window_idx] = copy(redunhaps_bitvec1) # avoids false sharing
            redundant_haplotypes[k].strand2[window_idx] = copy(redunhaps_bitvec2)
        end
    end

    return nothing
end

"""
    haplopair(X, H)

Calculate the best pair of haplotypes in `H` for each individual in `X`. Missing data in `X` 
does not have missing data. Missing data is initialized as 2x alternate allele freq.

# Input
* `X`: `p x n` genotype matrix possibly with missings. Each column is an individual.
* `H`: `p * d` haplotype matrix. Each column is a haplotype.

# Output
* `happair`: optimal haplotype pairs. `X[:, k] ≈ H[:, happair[1][k]] + H[:, happair[2][k]]`.
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
"""
function haplopair!(
    X::AbstractMatrix, # p × n
    H::AbstractMatrix; # p × d
    # preallocated vectors
    happair1::AbstractVector = ones(Int, size(X, 2)),      # length n 
    happair2::AbstractVector = ones(Int, size(X, 2)),      # length n
    hapscore::AbstractVector = zeros(Float32, size(X, 2)), # length n
    # preallocated matrices
    M     :: Matrix{Float32} = zeros(Float32, size(H, 2), size(H, 2)), # cannot be preallocated until Julia 2.0
    Xwork :: AbstractMatrix{Float32} = zeros(Float32, size(X, 1), size(X, 2)), # p × n
    Hwork :: AbstractMatrix{Float32} = convert(Matrix{Float32}, H),            # p × d
    N     :: AbstractMatrix{Float32} = zeros(Float32, size(X, 2), size(H, 2)), # n × d
    # Hwork :: ElasticArray{Float32} = convert(ElasticArrays{Float32}, H),            # p × d
    # N     :: ElasticArray{Float32} = ElasticArrays{Float32}(undef, size(X, 2), size(H, 2)), # n × d
    )
    p, n  = size(X)
    d     = size(H, 2)

    # reallocate matrices for last window (TODO: Hwork)
    if size(Xwork, 1) != p
        Xwork = zeros(Float32, p, n)
        # Hwork = ElasticArray{Float32}(undef, p, d)
    end

    # resize N
    # ElasticArrays.resize!(N, n, d)
    # ElasticArrays.resize!(Hwork, p, d)
    # copyto!(Hwork, H)

    # initializes missing
    initXfloat!(Xwork, X)

    t2, t3 = haplopair!(Xwork, Hwork, M, N, happair1, happair2, hapscore)
    t1 = t4 = 0 # no time spent on haplotype thinning or rescreening

    return t1, t2, t3, t4
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
    X::AbstractMatrix{Float32},
    H::AbstractMatrix{Float32},
    M::AbstractMatrix{Float32},
    N::AbstractMatrix{Float32},
    happair1::AbstractVector{Int32},
    happair2::AbstractVector{Int32},
    hapscore::AbstractVector{Float32}
    )

    p, n, d = size(X, 1), size(X, 2), size(H, 2)

    # assemble M (upper triangular only)
    t2 = @elapsed begin 
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
    end

    # computational routine
    t3 = @elapsed haplopair!(happair1, happair2, hapscore, M, N)

    # supplement the constant terms in objective
    t3 += @elapsed begin @inbounds for j in 1:n
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
* `happair`: optimal haplotype pair for each individual.
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
    fill!(hapmin, Inf32)

    @inbounds for k in 1:d, j in 1:k
        Mjk = M[j, k]
        # loop over individuals
        @simd for i in 1:n
            score = Mjk - N[i, j] - N[i, k]

            # keep best happair (original code)
            if score < hapmin[i]
                hapmin[i], happair1[i], happair2[i] = score, j, k
            end

            # keep all happairs that are equally good
            # if score < hapmin[i]
            #     empty!(happairs[i])
            #     push!(happairs[i], (j, k))
            #     hapmin[i] = score
            # elseif score == hapmin[i]
            #     push!(happairs[i], (j, k))
            # end

            # keep happairs that within some range of best pair (but finding all of them requires a 2nd pass)
            # if score < hapmin[i]
            #     empty!(happairs[i])
            #     push!(happairs[i], (j, k))
            #     hapmin[i] = score
            # elseif score <= hapmin[i] + tol && length(happairs[i]) < 100
            #     push!(happairs[i], (j, k))
            # end

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

"""
    initXfloat!(Xfloat, X)

Initializes the matrix `Xfloat` where missing values of matrix `X` by `2 x` allele frequency
and nonmissing entries of `X` are converted to type `Float32` for subsequent BLAS routines. 

# Input
* `X` is a `p x n` genotype matrix. Each column is an individual.
* `Xfloat` is the `p x n` matrix of X where missing values are filled by 2x allele frequency. 
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

    any(isnan, Xfloat) && error("Xfloat contains NaN during initialization! Shouldn't happen!")
    any(isinf, Xfloat) && error("Xfloat contains Inf during initialization! Shouldn't happen!")
    any(ismissing, Xfloat) && error("Xfloat contains Missing during initialization! Shouldn't happen!")

    return nothing
end

"""
    chunks(people, haplotypes)

Determines how many windows per chunk will be processed at once based on 
estimated memory. Total memory usage will be roughly 80% of total RAM.

# Inputs
- `d`: average number of unique haplotypes per window
- `td`: total number of haplotypes
- `p`: number of typed SNPs per window
- `n`: number of samples
- `threads`: number of threads (this affects `M` and `N` only)
- `Xbytes`: number of bytes to store genotype matrix
- `Hbytes`: number of bytes to store compressed haplotypes

# Output
- Maximum windows per chunk

# Memory intensive items:
- `M`: requires `d × d × 4` bytes per thread
- `N`: requires `d × p × 4` bytes per thread
- `redundant_haplotypes`: requires `windows × 2td × n` bits where `windows` is number of windows per chunk
"""
function nchunks(
    d::Int, 
    td::Int, 
    p::Int,
    n::Int,
    threads::Int = Threads.nthreads(),
    Xbytes::Int = 0,
    compressed_Hunique::Union{Nothing, CompressedHaplotypes} = nothing
    )
    # system info
    system_memory_gb = Sys.total_memory() / 2^30
    system_memory_bits = 8000000000 * system_memory_gb
    usable_bits = round(Int, system_memory_bits * 0.8) # use 80% of total memory

    # estimate memory usage per window
    Mbits_per_win = 32d * d * threads
    Nbits_per_win = 32d * p * threads
    Rbits_per_win = 2 * td * n

    # calculate X and H's memory requirement in bits
    Xbits = 4Xbytes
    Hbits = 0
    if !isnothing(compressed_Hunique)
        # avoid computing sizes for vector of strings because they are slow
        Hbits += Base.summarysize(compressed_Hunique.CW)
        Hbits += Base.summarysize(compressed_Hunique.CW_typed)
        Hbits += Base.summarysize(compressed_Hunique.start)
        Hbits += Base.summarysize(compressed_Hunique.pos)
        Hbits += Base.summarysize(compressed_Hunique.altfreq)
    end

    return round(Int, (usable_bits - Hbits - Xbits - Nbits_per_win - Mbits_per_win) / Rbits_per_win)
end
