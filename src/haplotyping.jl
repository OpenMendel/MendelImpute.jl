using NullableArrays

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
#    for i in 1:n
#        j, k = happair[1][i], happair[2][i]
#        hapmin[i] = M[j, k] - N[i, j] - N[i, k]
#    end
    fill!(hapmin, typemax(eltype(hapmin)))
    # TODO: parallel computing
    @inbounds for k in 1:d, j in 1:k
        # loop over individuals
        @simd for i in 1:n
            score = M[j, k] - N[i, j] - N[i, k]
            if score < hapmin[i]
                hapmin[i], happair[1][i], happair[2][i] = score, j, k
            end
        end
    end
    return nothing

end

"""
    haplopair(X, H)

Calculate the best pair of haplotypes in `H` for each individual in `X`.

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

    n, p     = size(X)
    d        = size(H, 1)
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
* `X`: `p x n` genotype matrix. Each row is an individual.
* `H`: `p x d` haplotype matrix. Each row is a haplotype.
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
    # assemble M
    At_mul_B!(M, H, H)
    for j in 1:d, i in 1:(j - 1) # off-diagonal
        M[i, j] = 2M[i, j] + M[i, i] + M[j, j]
    end
    for j in 1:d # diagonal
        M[j, j] *= 4
    end
    # assemble N
    At_mul_B!(N, X, H)
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
    fillmissing!(X, H, haplopair)

Fill in missing genotypes in `X` according to haplotypes. Non-missing genotypes
remain same.

# Input
* `X`: `p x n` genotype matrix. Each column is an individual.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `happair`: pair of haplotypes. `X[:, k] = H[:, happair[1][k]] + H[:, happair[2][k]]`.

# Output
* `discrepancy`: sum of squared errors between current values in missing genotypes
    and the imputed genotypes.
"""
function fillmissing!(
    X::NullableMatrix,
    H::AbstractMatrix,
    happair::Tuple{AbstractVector, AbstractVector}
    )

    discrepancy = zero(promote_type(eltype(X.values), eltype(H)))
    @inbounds for j in 1:size(X, 2), i in 1:size(X, 1)
        if X.isnull[i, j]
            tmp = H[i, happair[1][j]] + H[i, happair[2][j]]
            discrepancy   += abs2(X.values[i, j] - tmp)
            X.values[i, j] = tmp
        end
    end
    return discrepancy

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
    initmissing(X)

Initialize the missing values in a nullable matrix `X` by `2 x` allele frequency.
`X` is a `p x n` genotype matrix. Each column is an individual.
"""
function initmissing!(X::NullableMatrix)
    T = eltype(X.values)
    p, n = size(X)
    for i in 1:p
        # allele frequency
        cnnz = 0
        csum = zero(T)
        for j in 1:n
            if !X.isnull[i, j]
                cnnz += 1
                csum += X.values[i, j]
            end
        end
        # set missing values to 2freq
        imp = csum / cnnz
        for j in 1:n
            if X.isnull[i, j]
                X.values[i, j] = imp
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
* `X`: `p x n` nullable matrix. Each column is genotypes of an individual.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `M`: overwritten by `M[i, j] = 2dot(H[:, i], H[:, j]) + sumabs2(H[:, i]) +
    sumabs2(H[:, j])`.
* `N`: overwritten by `n x d` matrix `2X'H`.
* `happair`: optimal haplotype pair. `X[:, k] ≈ H[:, happair[k, 1]] + H[:, happair[k, 2]]`.
* `hapscore`: haplotyping score. 0 means best. Larger means worse.
* `maxiters`: number of MM iterations. Defaultis 1.
* `tolfun`: convergence tolerance of MM iterations. Default is 1e-3.
"""
function haploimpute!(
    X::NullableMatrix,
    H::AbstractMatrix,
    M::AbstractMatrix,
    N::AbstractMatrix,
    happair::Tuple{AbstractVector, AbstractVector},
    hapscore::AbstractVector,
    maxiters::Int  = 1,
    tolfun::Number = 1e-3
    )

    obj = typemax(eltype(hapscore))
    initmissing!(X)
    for iter in 1:maxiters
        # haplotyping
        haplopair!(X.values, H, M, N, happair, hapscore)
        # impute missing entries according to current haplotypes
        discrepancy = fillmissing!(X, H, happair)
        #println("discrepancy = $discrepancy")
        # convergence criterion
        objold = obj
        obj = sum(hapscore) - discrepancy
        #println("iter = $iter, obj = $obj")
        if abs(obj - objold) < tolfun * (objold + 1)
            break
        end
    end
    return nothing

end

"""
    phase(X, H, width=400, verbose=true)

Phasing (haplotying) of genotype matrix `X` from a pool of haplotypes `H`
by sliding windows.

# Input
* `X`: `p x n` nullable matrix. Each column is genotypes of an individual.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `width`: width of the sliding window.
* `verbose`: display algorithmic information.
"""
function phase(
    X::NullableMatrix{T},
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
    Xwork = X[1:3width, :] # NullableMatrix
    Hwork = view(H, 1:3width, :)
    happair_prev = deepcopy(happair)

    # number of windows
    windows = floor(Int, snps / width)

    # phase and impute window 1
    if verbose; println("Imputing SNPs 1:$width"); end
    haploimpute!(Xwork, Hwork, M, N, happair, hapscore)
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
        copy!(Xwork, view(X, ((w - 2) * width + 1):((w + 1) * width), :))
        # overwrite first 1/3 according to phased haplotypes
        #Hw1 = view(Hwork, 1:width, :)
        #fillgeno!(Xw1, Hw1, happair)
        #fill!(Xwb1, false) # TODO DO THIS OR NOT
        # phase current window
        copy!(happair_prev[1], happair[1])
        copy!(happair_prev[2], happair[2])
        haploimpute!(Xwork, Hwork, M, N, happair, hapscore)
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

    # phase last window
    if verbose
        println("Imputing SNPs $((windows - 1) * width + 1):$snps")
    end
    Xwork = X[((windows - 2) * width + 1):snps, :]
    Hwork = view(H, ((windows - 2) * width + 1):snps, :)
    #Hw1   = view(Hwork, 1:width, :)
    #Xw1   = view(Xwork.values, 1:width, :)
    #fillgeno!(Xw1, Hw1, happair)
    copy!(happair_prev[1], happair[1])
    copy!(happair_prev[2], happair[2])
    haploimpute!(Xwork, Hwork, M, N, happair, hapscore)
    # find optimal break points and record info to phase
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
    # i | j
    # k | l
    if i == k && j == l
        return (i, j), (-1, -1)
    end
    # i | j
    # l | k
    if i == l && j == k
        return (i, j), (-1, -1)
    end

    tol = 5
    # only one strand matches
    # i | j
    # k | l
    if i == k && j ≠ l
        if columndiff(H, j, l) < tol
            return (i, j), (-1, -1)
        else
            breakpt, errors = search_breakpoint(X, H, i, (j, l))
            return (k, l), (-1, breakpt)
        end
    # i | j
    # l | k
    elseif i == l && j ≠ k
        if columndiff(H, j, k) < tol
            return (i, j), (-1, -1)
        else
            breakpt, errors = search_breakpoint(X, H, i, (j, k))
            return (i, k), (-1, breakpt)
        end
    # i | j
    # l | k
    elseif j == k && i ≠ l
        if columndiff(H, i, l) < tol
            return (i, j), (-1, -1)
        else
            breakpt, errors = search_breakpoint(X, H, j, (i, l))
            return (l, j), (breakpt, -1)
        end
    # i | j
    # k | l
    elseif j == l && i ≠ k
        if columndiff(H, i, k) < tol
            return (i, j), (-1, -1)
        else
            breakpt, errors = search_breakpoint(X, H, j, (i, k))
            return (k, j), (breakpt, -1)
        end
    end

    # no strand matches
    # i | j
    # k | l
    if columndiff(H, i, k) < tol # i == k
        if columndiff(H, j, l) < tol # j == l
            return (i, j), (-1, -1)
        else # j ≠ l
            breakpt, errors = search_breakpoint(X, H, i, (j, l))
            return (i, l), (-1, breakpt)
        end
    else # i ≠ k
        if columndiff(H, j, l) < tol # j == l
            breakpt, errors = search_breakpoint(X, H, j, (i, k))
            return (k, j), (breakpt, -1)
        else # j ≠ l
            # breakpts1, errors1 = search_breakpoint(X, H, (i, k), (j, l))
            return (k, l), (0, 0)
        end
    end
    # i | j
    # l | k
    if columndiff(H, i, l) < tol # i == l
        if columndiff(H, j, k) < tol # j == k
            return (i, j), (-1, -1)
        else # j ≠ k
            breakpt, errors = search_breakpoint(X, H, i, (j, k))
            return (i, k), (-1, breakpt)
        end
    else # i ≠ l
        if columndiff(H, j, k) < tol # j == k
            breakpt, errors = search_breakpoint(X, H, j, (i, l))
            return (l, -1), (breakpt, -1)
        else # j ≠ k
            # breakpts2, errors2 = search_breakpoint(X, H, (i, l), (j, k))
            return (l, k), (0, 0)
        end
    end
    # # choose the best one
    # if errors1 < errors2
    #     return (k, l), breakpts2
    # else
    #     return (l, k), breakpts2
    # end
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
        if !isnull(X[pos])
            errors += unsafe_get(X[pos]) ≠ H[pos, s1] + H[pos, s2[2]]
        end
    end
    bkpt_optim, err_optim = 0, errors
    # quick return if perfect match
    err_optim == 0 && return 0, 0
    # extend haplotype H[:, s2[1]] position by position
    @inbounds for bkpt in 1:n
        if !isnull(X[bkpt]) && H[bkpt, s2[1]] ≠ H[bkpt, s2[2]]
            errors -= unsafe_get(X[bkpt]) ≠ H[bkpt, s1] + H[bkpt, s2[2]]
            errors += unsafe_get(X[bkpt]) ≠ H[bkpt, s1] + H[bkpt, s2[1]]
            if errors < err_optim
                bkpt_optim, err_optim = bkpt, errors
                # quick return if perfect match
                err_optim == 0 && return bkpt_optim, err_optim
            end
        end
    end
    return bkpt_optim, err_optim

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
            if !isnull(X[pos])
                errors += unsafe_get(X[pos]) ≠ H[pos, s1[1]] + H[pos, s2[2]]
            end
        end
        for pos in (bkpt1 + 1):length(X)
            if !isnull(X[pos])
                errors += unsafe_get(X[pos]) ≠ H[pos, s1[2]] + H[pos, s2[2]]
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
            if !isnull(X[bkpt2])
                errors -= unsafe_get(X[bkpt2]) ≠ H[bkpt2, s1[1]] + H[bkpt2, s2[2]]
                errors += unsafe_get(X[bkpt2]) ≠ H[bkpt2, s1[1]] + H[bkpt2, s2[1]]
                if errors < err_optim
                    err_optim = errors
                    bkpts_optim = (bkpt1, bkpt2)
                end
            end
        end
        for bkpt2 in (bkpt1 + 1):length(X)
            if !isnull(X[bkpt2])
                errors -= unsafe_get(X[bkpt2]) ≠ H[bkpt2, s1[2]] + H[bkpt2, s2[2]]
                errors += unsafe_get(X[bkpt2]) ≠ H[bkpt2, s1[2]] + H[bkpt2, s2[1]]
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

"""
    columndiff(M, j, k, [p=2])

L1 norm of the difference between the `j`-th and `k`-th columns of a matrix `M`.
"""
function columndiff(M::AbstractMatrix, j::Integer, k::Integer)
    d = zero(eltype(M))
    @inbounds @simd for i in 1:size(M, 1)
        d += abs(M[i, j] - M[i, k])
    end
    d
end
