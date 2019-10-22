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
    trans::Char
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

    # record unique haplotypes and mappings window by window
    for w in 1:(windows-1)
        H_cur_window = view(H, ((w - 1) * width + 1):(w * width), :)
        groupslices!(hapset.hapmap[w], H_cur_window, dim)
        hapset.uniqueindex[w] = unique(hapset.hapmap[w])
    end

    # find unique haplotype in last window
    H_last_window = view(H, ((windows - 1) * width + 1):p, :)
    groupslices!(hapset.hapmap[end], H_last_window, dim)
    hapset.uniqueindex[end] = unique(hapset.hapmap[end])

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
    X::AbstractMatrix{Union{Missing, T}},
    H::AbstractMatrix{T};
    width::Int    = 128,
    verbose::Bool = true
    ) where T <: Real

    # define some constants
    snps, people = size(X)
    haplotypes = size(H, 2)
    windows = floor(Int, snps / width)

    # get unique haplotype indices and maps for each window
    Hunique  = unique_haplotypes(H, width, 'T')
    num_uniq = length(Hunique.uniqueindex[1])

    # Initialize data structure for redundant haplotypes that matches the optimal one. 
    optimal_haplotypes = [OptimalHaplotypeSet(windows, haplotypes) for i in 1:people]

    # allocate working arrays
    happair     = ones(Int, people), ones(Int, people)
    hapscore    = zeros(T, people)
    Hwork       = ElasticArray{T}(H[1:width, Hunique.uniqueindex[1]])
    Xwork       = X[1:width, :]
    Xwork_float = zeros(T, size(Xwork))
    M           = zeros(T, num_uniq, num_uniq)
    N           = ElasticArray{T}(undef, people, num_uniq)

    # In first window, calculate optimal haplotype pair among unique haplotypes
    haploimpute!(Xwork, Hwork, M, N, happair, hapscore, Xfloat=Xwork_float)

    # find all haplotypes matching the optimal haplotype pairs
    compute_redundant_haplotypes!(optimal_haplotypes, Hunique, happair, H, 1)

    #TODO: make this loop multithreaded 
    for w in 2:(windows-1)

        # sync Xwork and Hwork with original data
        cur_range = ((w - 1) * width + 1):(w * width)
        M = resize_and_sync!(Xwork, Hwork, Hunique.uniqueindex[w], cur_range, X, H, M, N)

        # Calculate optimal haplotype pair among unique haplotypes
        haploimpute!(Xwork, Hwork, M, N, happair, hapscore, Xfloat=Xwork_float)

        # find all haplotypes matching the optimal haplotype pairs
        compute_redundant_haplotypes!(optimal_haplotypes, Hunique, happair, H, w)
    end

    # last window
    last_range = ((windows - 1) * width + 1):snps
    if mod(length(last_range), width) == 0
        #resize the typical way if the last window has the same width as previous windows
        M = resize_and_sync!(Xwork, Hwork, Hunique.uniqueindex[end], last_range, X, H, M, N)
        haploimpute!(Xwork, Hwork, M, N, happair, hapscore)
        compute_redundant_haplotypes!(optimal_haplotypes, Hunique, happair, H, windows)
    else
        #reallocate everything 
        num_uniq    = length(Hunique.uniqueindex[end])
        Hwork       = H[last_range, Hunique.uniqueindex[end]]
        Xwork       = X[last_range, :]
        M           = zeros(T, num_uniq, num_uniq)
        N           = zeros(T, people, num_uniq)
        haploimpute!(Xwork, Hwork, M, N, happair, hapscore)
        compute_redundant_haplotypes!(optimal_haplotypes, Hunique, happair, H, windows)
    end

    return optimal_haplotypes
end

"""
Helper function for recording optimal-redundant haplotypes for each window.
"""
function compute_redundant_haplotypes!(
    optimal_haplotypes::Vector{OptimalHaplotypeSet}, 
    Hunique::UniqueHaplotypeMaps, 
    happair::Tuple{AbstractVector, AbstractVector}, 
    H::AbstractMatrix,
    window::Int,
    )

    people = length(optimal_haplotypes)

    # loop through all people
    @inbounds for k in 1:people
        (Hwork_i, Hwork_j) = (happair[1][k], happair[2][k])
        # println("person $k's optimal haplotype pairs are: $((Hwork_i, Hwork_j))")

        (H_i, H_j) = (Hunique.uniqueindex[window][Hwork_i], Hunique.uniqueindex[window][Hwork_j])
        # println("person $k's optimal haplotype pairs are located at columns $H_i and $H_j in H")

        # loop through all haplotypes and find ones that match either of the optimal haplotypes 
        for jj in 1:size(H, 2)
            Hunique.hapmap[window][jj] == H_i && (optimal_haplotypes[k].strand1[window][jj] = true)
            Hunique.hapmap[window][jj] == H_j && (optimal_haplotypes[k].strand2[window][jj] = true)
        end
    end

    return nothing
end

"""
    resize_and_sync!(X, H, M, N, Xwork, Hwork, Hnext, window)

Up/downsizes the dimension of `Hwork`, `M`, and `N` and copies relevant information into `Xwork` and `Hwork`. 

# Inputs
* `Xwork`: Worker matrix storing X[window, :]. 
* `Hwork`: Haplotype matrix in the current window containing only unique haplotypes. Must add/subtract columns. 
* `Hnext`: The unique haplotype indices of the next haplotype window. 
* `window`: Indices of current window. 
* `X`: Full genotype matrix. Each column is a person's haplotype
* `H`: Full haplotype reference panel. Each column is a haplotype
* `M`: Square matrix used in the computational routine. Must be resized in both dimension. 
* `N`: Matrix used in the computational routine. Must add/subtract columns. 

TODO: check how ReshapedArray cause type instability and if it is significant overhead
"""
function resize_and_sync!(
    Xwork::AbstractMatrix,
    Hwork::ElasticArray,
    Hnext::Vector{Int},
    window::UnitRange{Int},
    X::AbstractMatrix,
    H::AbstractMatrix,
    M::AbstractMatrix,
    N::ElasticArray,
    )

    pp, dd = size(Hwork)
    next_d = length(Hnext)

    # resize working arrays
    if dd != next_d
        resize!(Hwork, pp        , next_d)
        resize!(N    , size(N, 1), next_d)
        Mvec = vec(M)
        resize!(Mvec, next_d^2)
        Mnew = Base.ReshapedArray(Mvec, (next_d, next_d), ()) # actually resize! makes a copy internally!
        # Mnew = zeros(eltype(M), next_d, next_d)               # always reallocate entire M
        # Mnew = (next_d < dd ? Base.ReshapedArray(vec(M), (next_d, next_d), ()) : 
        #                       zeros(eltype(M), next_d, next_d))
    else
        Mnew = M
    end

    # sync Xwork and Hwork with original data
    copyto!(Xwork, view(X, window, :))
    copyto!(Hwork, view(H, window, Hnext))

    return Mnew
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
    strand1::BitVector,
    strand2::Tuple{BitVector, BitVector}
    )

    n = length(X)

    # all haplotypes in BitSet are equivalent in current window, so get one as representative
    s1  = findfirst(strand1)
    s21 = findfirst(strand2[1])
    s22 = findfirst(strand2[2])

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

"""
    groupslices(A, dim)

Returns a vector of integers where each integer element of the returned vector
is a group number corresponding to the unique slices along dimension `dim` as
returned from `unique(A, dim)`, where `A` can be a multidimensional array.

# Example usage:
If `C = unique(A, dim)`, `ic = groupslices(A, dim)`, and
`ndims(A) == ndims(C) == 3`, then:
```
if dim == 1
   all(A .== C[ic,:,:])
elseif dim == 2
   all(A .== C[:,ic,:])
elseif dim == 3
   all(A .== C[:,:,ic])
end
```

Function from: https://github.com/mcabbott/GroupSlices.jl/blob/master/src/GroupSlices.jl
Can delete this function when this issue gets resolved: https://github.com/JuliaLang/julia/issues/1845 
"""
@generated function groupslices(A::AbstractArray{T,N}, dim::Int) where {T,N}
    quote
        if !(1 <= dim <= $N)
            ArgumentError("Input argument dim must be 1 <= dim <= $N, but is currently $dim")
        end
        hashes = zeros(UInt, size(A, dim))

        # Compute hash for each row
        k = 0
        @nloops $N i A d->(if d == dim; k = i_d; end) begin
            @inbounds hashes[k] = hash(hashes[k], hash((@nref $N A i)))
        end

        # Collect index of first row for each hash
        uniquerow = Vector{Int}(undef, size(A, dim))
        firstrow = Dict{Prehashed,Int}()
        for k = 1:size(A, dim)
            uniquerow[k] = get!(firstrow, Prehashed(hashes[k]), k)
        end
        uniquerows = collect(values(firstrow))

        # Check for collisions
        collided = falses(size(A, dim))
        @inbounds begin
            @nloops $N i A d->(if d == dim
                k = i_d
                j_d = uniquerow[k]
            else
                j_d = i_d
            end) begin
                if (@nref $N A j) != (@nref $N A i)
                    collided[k] = true
                end
            end
        end

        if any(collided)
            nowcollided = BitArray(size(A, dim))
            while any(collided)
                # Collect index of first row for each collided hash
                empty!(firstrow)
                for j = 1:size(A, dim)
                    collided[j] || continue
                    uniquerow[j] = get!(firstrow, Prehashed(hashes[j]), j)
                end
                for v in values(firstrow)
                    push!(uniquerows, v)
                end

                # Check for collisions
                fill!(nowcollided, false)
                @nloops $N i A d->begin
                    if d == dim
                        k = i_d
                        j_d = uniquerow[k]
                        (!collided[k] || j_d == k) && continue
                    else
                        j_d = i_d
                    end
                end begin
                    if (@nref $N A j) != (@nref $N A i)
                        nowcollided[k] = true
                    end
                end
                (collided, nowcollided) = (nowcollided, collided)
            end
        end
        ie = unique(uniquerow)
        ic_dict = Dict{Int,Int}()
        for k = 1:length(ie)
            ic_dict[ie[k]] = k
        end

        ic = similar(uniquerow)
        for k = 1:length(ic)
            ic[k] = ie[ic_dict[uniquerow[k]]]
        end
        return ic
    end
end

@generated function groupslices!(ic::AbstractArray, A::AbstractArray{T,N}, dim::Int) where {T,N}
    quote
        if !(1 <= dim <= $N)
            ArgumentError("Input argument dim must be 1 <= dim <= $N, but is currently $dim")
        end
        hashes = zeros(UInt, size(A, dim))

        # Compute hash for each row
        k = 0
        @nloops $N i A d->(if d == dim; k = i_d; end) begin
            @inbounds hashes[k] = hash(hashes[k], hash((@nref $N A i)))
        end

        # Collect index of first row for each hash
        uniquerow = Vector{Int}(undef, size(A, dim))
        firstrow = Dict{Prehashed,Int}()
        for k = 1:size(A, dim)
            uniquerow[k] = get!(firstrow, Prehashed(hashes[k]), k)
        end
        uniquerows = collect(values(firstrow))

        # Check for collisions
        collided = falses(size(A, dim))
        @inbounds begin
            @nloops $N i A d->(if d == dim
                k = i_d
                j_d = uniquerow[k]
            else
                j_d = i_d
            end) begin
                if (@nref $N A j) != (@nref $N A i)
                    collided[k] = true
                end
            end
        end

        if any(collided)
            nowcollided = BitArray(size(A, dim))
            while any(collided)
                # Collect index of first row for each collided hash
                empty!(firstrow)
                for j = 1:size(A, dim)
                    collided[j] || continue
                    uniquerow[j] = get!(firstrow, Prehashed(hashes[j]), j)
                end
                for v in values(firstrow)
                    push!(uniquerows, v)
                end

                # Check for collisions
                fill!(nowcollided, false)
                @nloops $N i A d->begin
                    if d == dim
                        k = i_d
                        j_d = uniquerow[k]
                        (!collided[k] || j_d == k) && continue
                    else
                        j_d = i_d
                    end
                end begin
                    if (@nref $N A j) != (@nref $N A i)
                        nowcollided[k] = true
                    end
                end
                (collided, nowcollided) = (nowcollided, collided)
            end
        end
        ie = unique(uniquerow)
        ic_dict = Dict{Int,Int}()
        for k = 1:length(ie)
            ic_dict[ie[k]] = k
        end

        for k = 1:length(ic)
            ic[k] = ie[ic_dict[uniquerow[k]]]
        end
    end
end

# """
#     unique_haplotypes(H, width)

# Finds the unique haplotypes determined by the reference haplotypes stored 
# in the columns of H. 

# Note: This seems to be ~3x slower than using groupslices. 

# # Input
# * `H`: an `p x d` reference panel of haplotypes. 
# * `width`: The window width 
# """
# function unique_haplotypes(
#     H::AbstractMatrix, 
#     width::Int,
#     )

#     p, d    = size(H)
#     windows = ceil(Int, p / width)
#     fast_data_type = Dict(8=>UInt8, 16=>UInt16, 32=>UInt32, 64=>UInt64, 128=>UInt128)

#     if eltype(H) == Bool && haskey(fast_data_type, width)
#         return fast_elimination(H, windows, width, H[1:width, :], fast_data_type)
#     else
#         ??
#     end
# end

# """
#     fast_elimination!(unique_hap, H, windows, width)

# Computes the columns of `H` that are unique in each window and stores non-unique mappings. 

# # Input
# * `H`: an `p x d` reference panel of haplotypes. 
# * `windows`: total number of windows
# * `width`: the width of a window (should be 8, 16, 32, 64, or 128)
# * `storage`: an `width x d` Bitmatrix 
# * `fast_data_type`: the data types that can use fast_elimination

# # Output
# * `unique_hap` that stores the correct unique haplotypes and mappings for non-unique haplotypes
# """
# function fast_elimination(
#     H::BitMatrix, 
#     windows::Int64, 
#     width::Int64,
#     storage::BitMatrix = H[1:width, :],
#     fast_data_type::Dict = Dict(8=>UInt8, 16=>UInt16, 32=>UInt32, 64=>UInt64, 128=>UInt128)
#     )

#     unique_hap = UniqueHaplotypes(windows, size(H, 2))

#     # reinterpret each haplotype as an integer
#     HR = reinterpret(fast_data_type[width], storage.chunks) 

#     # record unique haplotypes and non-unique mappings in first window
#     unique_index!(unique_hap.unique_index[1], unique_hap.redundant_map[1], HR)

#     # loop through windows
#     for w in 2:windows-1
#         copyto!(storage, @view(H[((w - 1) * width + 1):(w * width), :]))
#         HR = reinterpret(fast_data_type[width], storage.chunks) 
#         unique_index!(unique_hap.unique_index[w], unique_hap.redundant_map[w], HR)
#     end

#     # TODO: last window may have length ∉ fast_data_type
#     return nothing
# end

# # helper function for fast_elimination!
# function unique_index!(u::BitVector, d::Dict{Int64, Int64}, v::AbstractVector)
#     seen = Set{eltype(v)}()

#     @inbounds for i in 1:length(v)
#         if v[i] ∈ seen
#             u[i] = false
#             d[i] = findfirst(isequal(v[i]), v)
#         else
#             push!(seen, v[i])
#         end
#     end
# end
