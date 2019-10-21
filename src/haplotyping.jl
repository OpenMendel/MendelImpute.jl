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
    happair  = ones(Int, n), ones(Int, n)
    hapscore = zeros(eltype(N), n)
    haplopair!(X, H, M, N, happair, hapscore)

    return happair, hapscore
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
    happair::Tuple{AbstractVector, AbstractVector},
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
    haplopair!(happair, hapscore, M, N)

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

Calculate the best pair of haplotypes in `H` for each individual in `X` using
sufficient statistics `M` and `N`.

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
    happair::Tuple{AbstractVector, AbstractVector},
    hapmin::Vector,
    M::AbstractMatrix,
    N::AbstractMatrix
    )

    n, d = size(N)
    fill!(hapmin, typemax(eltype(hapmin)))

    @inbounds for k in 1:d, j in 1:k
        # loop over individuals
        for i in 1:n
            score = M[j, k] - N[i, j] - N[i, k]
            if score < hapmin[i]
                hapmin[i], happair[1][i], happair[2][i] = score, j, k
            end
        end
    end

    return nothing
end

"""
    fillmissing!(X, H, haplopair)

Fill in missing genotypes in `X` according to haplotypes. Non-missing genotypes
remain same.

# Input
* `X`: `p x n` genotype matrix. Each column is an individual.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `happair`: pair of haplotypes. `X[:, k] = H[:, happair[1][k]] + H[:, happair[2][k]]`.
"""
function fillmissing!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    happair::Tuple{AbstractVector, AbstractVector},
    )

    p, n = size(X)

    @inbounds for j in 1:n, i in 1:p
        if ismissing(X[i, j])
            X[i, j] = H[i, happair[1][j]] + H[i, happair[2][j]]
        end
    end

    return nothing
end

"""
    fillgeno!(X, H, happair)

Fill in genotypes according to haplotypes. Both missing and non-missing
genotypes may be changed.

# Input
* `X`: `p x n` genotype matrix. Each column is an individual.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `happair`: pair of haplotypes. `X[:, k] = H[:, happair[1][k]] + H[:, happair[2][k]]`.
"""
function fillgeno!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    happair::Tuple{AbstractVector, AbstractVector}
    )

    @inbounds for j in 1:size(X, 2), i in 1:size(X, 1)
        X[i, j] = H[i, happair[1][j]] + H[i, happair[2][j]]
    end
    return nothing

end

"""
    initmissing(X, Xwork)

Initializes the matrix `Xfloat` where missing values of matrix `X` by `2 x` allele frequency.

# Input
* `X` is a `p x n` genotype matrix. Each column is an individual.
* `Xfloat` is the `p x n` matrix of X where missing values are filled by 2x allele frequency. 
"""
function initmissing!(
    X::AbstractMatrix;
    Xfloat::AbstractMatrix = zeros(Float32, size(X))
    )
    
    T = eltype(X)
    p, n = size(X)

    for i in 1:p
        # allele frequency
        cnnz = 0
        csum = zero(T)
        for j in 1:n
            if !ismissing(X[i, j])
                cnnz += 1
                csum += X[i, j]
            end
        end
        # set missing values to 2freq
        imp = csum / cnnz
        for j in 1:n
            if ismissing(X[i, j]) 
                Xfloat[i, j] = imp
            else
                Xfloat[i, j] = X[i, j]
            end
        end
    end

    return nothing
end

"""
    haploimpute!(X, H, M, N, happair, hapscore, maxiters=1, tolfun=1e-3)

Haplotying of genotype matrix `X` from the pool of haplotypes `H` and impute
missing genotypes in `X` according to haplotypes.

# Input
* `X`: `p x n` matrix with missing values. Each column is genotypes of an individual.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `M`: overwritten by `M[i, j] = 2dot(H[:, i], H[:, j]) + sumabs2(H[:, i]) +
    sumabs2(H[:, j])`.
* `N`: overwritten by `n x d` matrix `2X'H`.
* `happair`: optimal haplotype pair. `X[:, k] ≈ H[:, happair[k, 1]] + H[:, happair[k, 2]]`.
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
* `Xfloat`: copy of `X` where missing values are filled with mean. This engages in linear algebra for computing `N`
* `maxiters`: number of MM iterations. Default is 1.
* `tolfun`: convergence tolerance of MM iterations. Default is 1e-3.
"""
function haploimpute!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    M::AbstractMatrix,
    N::AbstractMatrix,
    happair::Tuple{AbstractVector, AbstractVector},
    hapscore::AbstractVector;
    Xfloat::AbstractMatrix = zeros(eltype(M), size(X)),
    maxiters::Int  = 1,
    tolfun::Number = 1e-3
    )

    obj = typemax(eltype(hapscore))
    initmissing!(X, Xfloat=Xfloat) #Xfloat[i, j] = X[i, j] on observed entries

    # haplotyping
    haplopair!(Xfloat, H, M, N, happair, hapscore)

    # impute missing entries according to current haplotypes
    # fillmissing!(X, H, happair)

    return nothing
end

"""
    phase(X, H, width=400, verbose=true)

Phasing (haplotying) of genotype matrix `X` from a pool of haplotypes `H`
by sliding windows.

# Input
* `X`: `p x n` matrix with missing values. Each column is genotypes of an individual.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `width`: width of the sliding window.
* `verbose`: display algorithmic information.
"""
function phase(
    X::AbstractMatrix{Union{Missing, T}},
    H::AbstractMatrix{T},
    width::Int    = 400,
    verbose::Bool = true
    ) where T <: Real

    people, snps, haplotypes = size(X, 2), size(X, 1), size(H, 2)
    # allocate working arrays
    M        = zeros(T, haplotypes, haplotypes)
    N        = zeros(T,     people, haplotypes)
    happair  = ones(Int, people), ones(Int, people)
    hapscore = zeros(T, people)
    phase    = [HaplotypeMosaicPair(snps) for i in 1:people]

    # no need for sliding window
    if snps ≤ 3width
        haploimpute!(X, H, M, N, happair, hapscore)
        for i in 1:people
            push!(phase[i].strand1.start, 1)
            push!(phase[i].strand1.haplotypelabel, happair[1][i])
            push!(phase[i].strand2.start, 1)
            push!(phase[i].strand2.haplotypelabel, happair[2][i])
        end
        return phase
    end

    # allocate working arrays
    Xwork = X[1:3width, :]
    Xwork_float = zeros(T, size(Xwork))
    Hwork = view(H, 1:3width, :)
    # Hwork, (pp, dd) = unique_haplotypes(H, ((windows - 2) * width + 1):snps)
    happair_prev = deepcopy(happair)

    # number of windows
    windows = floor(Int, snps / width)

    # phase and impute window 1
    verbose && println("Imputing SNPs 1:$width")
    haploimpute!(Xwork, Hwork, M, N, happair, hapscore, Xfloat=Xwork_float)
    for i in 1:people
        push!(phase[i].strand1.start, 1)
        push!(phase[i].strand1.haplotypelabel, happair[1][i])
        push!(phase[i].strand2.start, 1)
        push!(phase[i].strand2.haplotypelabel, happair[2][i])
    end

    # first  1/3: ((w - 2) * width + 1):((w - 1) * width)
    # middle 1/3: ((w - 1) * width + 1):(      w * width)
    # last   1/3: (      w * width + 1):((w + 1) * width)
    for w in 2:(windows - 1)
        if verbose
            println("Imputing SNPs $((w - 1) * width + 1):$(w * width)")
        end

        # sync Xwork and Hwork with original data
        Hwork = view(H, ((w - 2) * width + 1):((w + 1) * width), :)
        # Hwork, (pp, dd) = unique_haplotypes(H, ((w - 2) * width + 1):((w + 1) * width))
        copyto!(Xwork, view(X, ((w - 2) * width + 1):((w + 1) * width), :))

        # phase current window
        copyto!(happair_prev[1], happair[1])
        copyto!(happair_prev[2], happair[2])
        haploimpute!(Xwork, Hwork, M, N, happair, hapscore, Xfloat=Xwork_float)
        # haploimpute!(Xwork, Hwork, view(M, 1:dd, 1:dd), view(N, :, 1:dd), happair, hapscore, Xfloat=Xwork_float)

        # find optimal break points and record info into phase
        Hw12 = view(Hwork, 1:2width, :)
        for i in 1:people
            Xi = view(Xwork, 1:2width, i)
            (happair[1][i], happair[2][i]), bkpts =
                continue_haplotype(Xi, Hw12,
                (happair_prev[1][i], happair_prev[2][i]),
                (     happair[1][i],      happair[2][i]))
            # strand 1
            if bkpts[1] > -1 && bkpts[1] < 2width
                push!(phase[i].strand1.start, (w - 2) * width + 1 + bkpts[1])
                push!(phase[i].strand1.haplotypelabel, happair[1][i])
            end
            # strand 2
            if bkpts[2] > -1 && bkpts[2] < 2width
                push!(phase[i].strand2.start, (w - 2) * width + 1 + bkpts[2])
                push!(phase[i].strand2.haplotypelabel, happair[2][i])
            end
            # # for debug
            if verbose == true && i == 1
                println("happair = ($(happair[1][i]), $(happair[2][i]))")
                println("bkpts = $bkpts")
            end
        end
    end

    # Hua's code without searching breakpoints
    # for w in 2:(windows-1)
    #     if verbose
    #         println("Imputing SNPs $((w - 1) * width + 1):$(w * width)")
    #     end

    #     # sync Xwork and Hwork with original data
    #     Hwork = view(H, ((w - 2) * width + 1):((w + 1) * width), :)
    #     # Hwork, (pp, dd) = unique_haplotypes(H, ((w - 2) * width + 1):((w + 1) * width))
    #     copyto!(Xwork, view(X, ((w - 2) * width + 1):((w + 1) * width), :))

    #     # phase current window
    #     copyto!(happair_prev[1], happair[1])
    #     copyto!(happair_prev[2], happair[2])
    #     haploimpute!(Xwork, Hwork, M, N, happair, hapscore, Xfloat=Xwork_float)
    #     # haploimpute!(Xwork, Hwork, view(M, 1:dd, 1:dd), view(N, :, 1:dd), happair, hapscore, Xfloat=Xwork_float)

    #     # record info into phase
    #     for i in 1:people
    #         push!(phase[i].strand1.start, (w - 2) * width + 1)
    #         push!(phase[i].strand1.haplotypelabel, happair[1][i])
    #         push!(phase[i].strand2.start, (w - 2) * width + 1)
    #         push!(phase[i].strand2.haplotypelabel, happair[2][i])
    #     end
    # end

    # phase last window
    if verbose
        println("Imputing SNPs $((windows - 1) * width + 1):$snps")
    end
    Xwork = X[((windows - 2) * width + 1):snps, :]
    Hwork = view(H, ((windows - 2) * width + 1):snps, :)
    # Hwork, (pp, dd) = unique_haplotypes(H, ((windows - 2) * width + 1):snps)
    copyto!(happair_prev[1], happair[1])
    copyto!(happair_prev[2], happair[2])
    # haploimpute!(Xwork, Hwork, view(M, 1:dd, 1:dd), view(N, :, dd), happair, hapscore)
    haploimpute!(Xwork, Hwork, M, N, happair, hapscore)

    # # find optimal break points and record info to phase
    for i in 1:people
        (happair[1][i], happair[2][i]), bkpts =
        continue_haplotype(Xwork[:, i], Hwork,
            (happair_prev[1][i], happair_prev[2][i]),
            (happair[1][i], happair[2][i]))
        # strand 1
        if bkpts[1] > -1 && bkpts[1] < 2width
            push!(phase[i].strand1.start, (windows - 2) * width + 1 + bkpts[1])
            push!(phase[i].strand1.haplotypelabel, happair[1][i])
        end
        # strand 2
        if bkpts[2] > -1 && bkpts[2] < 2width
            push!(phase[i].strand2.start, (windows - 2) * width + 1 + bkpts[2])
            push!(phase[i].strand2.haplotypelabel, happair[2][i])
        end
    end

    return phase
end

function phase2(
    X::AbstractMatrix{Union{Missing, T}},
    H::AbstractMatrix{T};
    width::Int    = 128,
    verbose::Bool = true
    ) where T <: Real

    # problem dimensions
    snps, people = size(X)

    # number of windows
    windows = floor(Int, snps / width)

    # get redundant haplotype sets. 
    hapset = compute_optimal_halotype_set(X, H, width=width, verbose=verbose)

    # allocate working arrays
    phase = [HaplotypeMosaicPair(snps) for i in 1:people]
    store = ([copy(hapset.strand1[1, i]) for i in 1:people], [copy(hapset.strand2[1, i]) for i in 1:people])
    window_span = (ones(Int, people), ones(Int, people))

    # TODO: parallel computing
    # TODO: replace `intersect` and `intersect!` with fast set intersection using bisection/seesaw search
    @inbounds for i in 1:people, w in 2:windows

        # decide how to concatenate next 2 windows to previous windows based on the larger intersection
        A = intersect(store[1][i], hapset.strand1[w, i])
        B = intersect(store[1][i], hapset.strand2[w, i])
        if length(A) >= length(B)
            # no need to cross over
            a = A 
        else
            # cross over
            a = B
            hapset.strand1.p[w, i], hapset.strand2.p[w, i] = hapset.strand2.p[w, i], hapset.strand1.p[w, i]
        end
        b = intersect(store[2][i], hapset.strand2[w, i])

        # strand 1
        if isempty(a)
            # delete all nonmatching haplotypes in previous windows
            for ww in (w - window_span[1][i]):(w - 1)
                hapset.strand1.p[ww, i] = copy(store[1][i]) 
            end

            # update counters and storage
            store[1][i] = copy(hapset.strand1[w, i])
            window_span[1][i] = 1
        else
            intersect!(store[1][i], hapset.strand1[w, i])
            window_span[1][i] += 1
        end

        # strand 2
        if isempty(b)
            # delete all nonmatching haplotypes in previous windows
            for ww in (w - window_span[2][i]):(w - 1)
                hapset.strand2.p[ww, i] = copy(store[2][i]) 
            end

            # update counters and storage
            store[2][i] = copy(hapset.strand2[w, i])
            window_span[2][i] = 1
        else
            intersect!(store[2][i], hapset.strand2[w, i])
            window_span[2][i] += 1
        end
    end

    # TODO: there's a bug in computing redundant haplotypes since last window never agrees with 2nd to last window
    # handle last few windows separately, since they may not hit the isempty command
    for i in 1:people
        for ww in (windows - window_span[1][i] + 1):windows
            hapset.strand1.p[ww, i] = copy(store[1][i]) 
        end

        for ww in (windows - window_span[2][i] + 1):windows
            hapset.strand2.p[ww, i] = copy(store[2][i]) 
        end
    end

    # phase window 1
    for i in 1:people
        push!(phase[i].strand1.start, 1)
        push!(phase[i].strand1.haplotypelabel, first(hapset.strand1[1, i]))
        push!(phase[i].strand2.start, 1)
        push!(phase[i].strand2.haplotypelabel, first(hapset.strand2[1, i]))
    end

    #phase window by window without checking breakpoints
    for i in 1:people, w in 2:windows
        hap1 = first(hapset.strand1[w, i])
        hap2 = first(hapset.strand2[w, i])

        # strand 1
        push!(phase[i].strand1.start, (w - 1) * width + 1)
        push!(phase[i].strand1.haplotypelabel, hap1)

        # strand 2
        push!(phase[i].strand2.start, (w - 1) * width + 1)
        push!(phase[i].strand2.haplotypelabel, hap2)
    end

    # find optimal break points and record info to phase. 
    # store = ([copy(hapset.strand1[1, i]) for i in 1:people], [copy(hapset.strand2[1, i]) for i in 1:people])
    # for i in 1:people, w in 2:windows
        
    #     a = intersect(store[1][i], hapset.strand1[w, i])
    #     b = intersect(store[2][i], hapset.strand2[w, i])

    #     if isempty(a)
    #         # search breakpoints
    #         Xi = view(X, ((w - 2) * width + 1):(w * width), i)
    #         Hi = view(H, ((w - 2) * width + 1):(w * width), :)
    #         prev_and_cur_haplotypes = (hapset.strand1[w - 1, i], hapset.strand1[w, i])
    #         bkpt, hap, err_optim = search_breakpoint(Xi, Hi, hapset.strand2[w, i], prev_and_cur_haplotypes)

    #         # record info into phase
    #         push!(phase[i].strand1.start, (w - 2) * width + 1 + bkpt)
    #         push!(phase[i].strand1.haplotypelabel, hap)

    #         # update storage
    #         store[1][i] = copy(hapset.strand1[w, i])
    #     end

    #     if isempty(b)
    #         # search breakpoints
    #         Xi = view(X, ((w - 2) * width + 1):(w * width), i)
    #         Hi = view(H, ((w - 2) * width + 1):(w * width), :)
    #         prev_and_cur_haplotypes = (hapset.strand2[w - 1, i], hapset.strand2[w, i])
    #         bkpt, hap, err_optim = search_breakpoint(Xi, Hi, hapset.strand1[w, i], prev_and_cur_haplotypes)

    #         # record info into phase
    #         push!(phase[i].strand2.start, (w - 2) * width + 1 + bkpt)
    #         push!(phase[i].strand2.haplotypelabel, hap)

    #         # update storage
    #         store[2][i] = copy(hapset.strand2[w, i])
    #     end
    # end

    # finally, fill in missing entries of X
    # impute!(X, H, phase)
    impute2!(X, H, phase)

    return hapset
    # return phase, hapset, bkpts
end

"""
    continue_haplotype(X, H, happair_prev, happair_next)

Find the optimal concatenated haplotypes from unordered haplotype pairs in two
consecutive windows.

# Input
* `X`: an `n` vector of genotypes with {0, 1, 2} entries
* `H`: an `n x d` reference panel of haplotypes with {0, 1} entries
* `happair_prev`: unordered haplotypes `(i, j)` in the first window
* `happair_next`: unordered haplotypes `(k, l)` in the second window

# Output
* `happair_next_optimal`: optimal ordered haplotypes in the second window
* `breakpt`: break points in the ordered haplotypes
"""
function continue_haplotype(
    X::AbstractVector,
    H::AbstractMatrix,
    happair_prev::Tuple{Int, Int},
    happair_next::Tuple{Int, Int}
    )

    i, j = happair_prev
    k, l = happair_next

    # both strands match
    if i == k && j == l
        return (k, l), (-1, -1)
    end

    if i == l && j == k
        return (l, k), (-1, -1)
    end

    # only one strand matches
    if i == k && j ≠ l
        breakpt, errors = search_breakpoint(X, H, i, (j, l))
        return (k, l), (-1, breakpt)
    elseif i == l && j ≠ k
        breakpt, errors = search_breakpoint(X, H, i, (j, k))
        return (l, k), (-1, breakpt)
    elseif j == k && i ≠ l
        breakpt, errors = search_breakpoint(X, H, j, (i, l))
        return (l, k), (breakpt, -1)
    elseif j == l && i ≠ k
        breakpt, errors = search_breakpoint(X, H, j, (i, k))
        return (k, l), (breakpt, -1)
    end

    return (k, l), (0, 0)

end

"""
    search_breakpoint(X, H, s1, s2)

Find the optimal break point between s2[1] and s2[2] in configuration
s1 | s2[1]
s1 | s2[2]
"""
function search_breakpoint(
    X::AbstractVector,
    H::AbstractMatrix,
    s1::Int,
    s2::Tuple{Int, Int}
    )

    n = length(X)
    # count number of errors if second haplotype is all from H[:, s2[2]]
    errors = 0
    for pos in 1:n
        if !ismissing(X[pos])
            errors += X[pos] ≠ H[pos, s1] + H[pos, s2[2]]
        end
    end
    bkpt_optim, err_optim = 0, errors

    # quick return if perfect match
    err_optim == 0 && return 0, 0

    # extend haplotype H[:, s2[1]] position by position
    @inbounds for bkpt in 1:n
        if !ismissing(X[bkpt]) && H[bkpt, s2[1]] ≠ H[bkpt, s2[2]]
            errors -= X[bkpt] ≠ H[bkpt, s1] + H[bkpt, s2[2]]
            errors += X[bkpt] ≠ H[bkpt, s1] + H[bkpt, s2[1]]
            if errors < err_optim
                bkpt_optim, err_optim = bkpt, errors
                # quick return if perfect match
                err_optim == 0 && return bkpt_optim, err_optim
            end
        end
    end

    return bkpt_optim, err_optim
end

function search_breakpoint(
    X::AbstractVector,
    H::AbstractMatrix,
    strand1::BitSet,
    strand2::Tuple{BitSet, BitSet}
    )

    n = length(X)

    # all haplotypes in BitSet are equivalent in current window, so get one as representative
    s1  = first(strand1)
    s21 = first(strand2[1])
    s22 = first(strand2[2])

    # count number of errors if second haplotype is all from H[:, s2[2]]
    errors = 0
    for pos in 1:n
        if !ismissing(X[pos])
            errors += X[pos] ≠ H[pos, s1] + H[pos, s22]
        end
    end
    bkpt_optim, err_optim = 0, errors

    # quick return if perfect match
    err_optim == 0 && return 0, s22, 0

    # extend haplotype H[:, s2[1]] position by position
    @inbounds for bkpt in 1:n
        if !ismissing(X[bkpt]) && H[bkpt, s21] ≠ H[bkpt, s22]
            errors -= X[bkpt] ≠ H[bkpt, s1] + H[bkpt, s22]
            errors += X[bkpt] ≠ H[bkpt, s1] + H[bkpt, s21]
            if errors < err_optim
                bkpt_optim, err_optim = bkpt, errors
                # quick return if perfect match
                err_optim == 0 && return bkpt_optim, s22, err_optim
            end
        end
    end

    return bkpt_optim, s22, err_optim
end

"""
    search_breakpoint(X, H, s1, s2)

Find the optimal break point between s2[1] and s2[2] in configuration
s1[1] | s2[1]
s1[2] | s2[2]
"""
function search_breakpoint(
    X::AbstractVector,
    H::AbstractMatrix,
    s1::Tuple{Int, Int},
    s2::Tuple{Int, Int}
    )

    err_optim   = typemax(Int)
    bkpts_optim = (0, 0)

    # search over all combintations of break points in two strands
    @inbounds for bkpt1 in 0:length(X)

        # count number of errors if second haplotype is all from H[:, s2[2]]
        errors = 0
        for pos in 1:bkpt1
            if !ismissing(X[pos])
                errors += X[pos] ≠ H[pos, s1[1]] + H[pos, s2[2]]
            end
        end
        for pos in (bkpt1 + 1):length(X)
            if !ismissing(X[pos])
                errors += X[pos] ≠ H[pos, s1[2]] + H[pos, s2[2]]
            end
        end
        if errors < err_optim
            err_optim = errors
            bkpts_optim = (bkpt1, 0)

            # quick return if perfect match
            err_optim == 0 && return bkpts_optim, err_optim
        end

        # extend haplotype H[:, s2[1]] position by position
        for bkpt2 in 1:bkpt1
            if !ismissing(X[bkpt2])
                errors -= X[bkpt2] ≠ H[bkpt2, s1[1]] + H[bkpt2, s2[2]]
                errors += X[bkpt2] ≠ H[bkpt2, s1[1]] + H[bkpt2, s2[1]]
                if errors < err_optim
                    err_optim = errors
                    bkpts_optim = (bkpt1, bkpt2)
                end
            end
        end
        for bkpt2 in (bkpt1 + 1):length(X)
            if !ismissing(X[bkpt2])
                errors -= X[bkpt2] ≠ H[bkpt2, s1[2]] + H[bkpt2, s2[2]]
                errors += X[bkpt2] ≠ H[bkpt2, s1[2]] + H[bkpt2, s2[1]]
                if errors < err_optim
                    err_optim = errors
                    bkpts_optim = (bkpt1, bkpt2)
                    # quick return if perfect match
                    err_optim == 0 && return bkpts_optim, err_optim
                end
            end
        end
    end

    return bkpts_optim, err_optim
end

function impute!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    phase::Vector{HaplotypeMosaicPair}
    )

    fill!(X, 0)
    # loop over individuals
    for i in 1:size(X, 2)
        for s in 1:(length(phase[i].strand1.start) - 1)
            idx = phase[i].strand1.start[s]:(phase[i].strand1.start[s + 1] - 1)
            X[idx, i] = H[idx, phase[i].strand1.haplotypelabel[s]]
        end
        idx = phase[i].strand1.start[end]:phase[i].strand1.length
        X[idx, i] = H[idx, phase[i].strand1.haplotypelabel[end]]
        for s in 1:(length(phase[i].strand2.start) - 1)
            idx = phase[i].strand2.start[s]:(phase[i].strand2.start[s + 1] - 1)
            X[idx, i] += H[idx, phase[i].strand2.haplotypelabel[s]]
        end
        idx = phase[i].strand2.start[end]:phase[i].strand2.length
        X[idx, i] += H[idx, phase[i].strand2.haplotypelabel[end]]
    end
end

function impute2!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    phase::Vector{HaplotypeMosaicPair}
    )

    p, n = size(X)

    @inbounds for person in 1:n, snp in 1:p
        if ismissing(X[snp, person])
            #find where snp is located in phase
            hap1_position = searchsortedlast(phase[person].strand1.start, snp)
            hap2_position = searchsortedlast(phase[person].strand2.start, snp)

            #find the correct haplotypes 
            hap1 = phase[person].strand1.haplotypelabel[hap1_position]
            hap2 = phase[person].strand2.haplotypelabel[hap2_position]

            # imputation step 
            X[snp, person] = H[snp, hap1] + H[snp, hap2]
        end
    end

    return nothing
end
