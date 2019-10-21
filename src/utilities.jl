
"""
    unique_haplotypes(H, window::UnitRange{Int})

Finds the unique haplotypes determined by the reference haplotypes stored 
in the columns of H. 

# Input
* `H`: an `p x d` reference panel of haplotypes within a genomic window. 
* `window`: a small window of `H` that is currently undergoing haplotyping.

# Output
* A `view` of `H` at the appropriate window with all redundant haplotypes eliminated
"""
function unique_haplotypes(
    H::AbstractMatrix, 
    window::UnitRange{Int}
    )

    lw = length(window)
    cur_chunk = view(H, window, :)

    if eltype(H) == Bool && lw in Set([8, 16, 32, 64, 128])
        unique_hap_index = unique_haplotype_idx(cur_chunk)
    else
        unique_hap_index = unique(groupslices(cur_chunk))
    end

    unique_hap_index = unique(groupslices(cur_chunk, 2))
    unique_hap = view(H, window, unique_hap_index)
    return unique_hap, size(unique_hap)
end

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

TODO: replace `groupslices!` with fast haplotype elimination strategy when width is a small multiple of 2
"""
function unique_haplotypes(
    H::AbstractMatrix,
    width::Int,
    trans::Char='N'
    )

    if trans == 'N'
        dim = 1
    elseif trans == 'T'
        dim = 2
    else
        error("trans can only be 'N' or 'T' but was $dim" )
    end

    p, d    = size(H)
    windows = ceil(Int, p / width)
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

# function unique_haplotypes(H::BitArray{2})
#     p, d = size(H) 

#     # reinterpret each haplotype as an integer
#     if p == 8 
#         HR = reinterpret(UInt8, H.chunks) 
#     elseif p == 16
#         HR = reinterpret(UInt16, H.chunks)
#     elseif p == 32
#         HR = reinterpret(UInt32, H.chunks)
#     elseif p == 64
#         HR = reinterpret(UInt64, H.chunks)
#     elseif p == 128
#         HR = reinterpret(UInt128, H.chunks)
#     else
#         return convert(Matrix{Float32}, unique(H, dims=1))
#     end
    
#     Hrank = denserank(HR) # map to unique integers with no gap
#     HU    = unique(HR)    # find unique integers
#     n     = length(HU)
#     Hrep  = zeros(Int, n) # representative haplotype for integer 

#     m = 0
#     for j = 1:d
#         if Hrep[Hrank[j]] == 0
#             Hrep[Hrank[j]] = j
#             m += 1
#             m == n && break
#         end
#     end

#     Hunique = convert(Matrix{Float32}, H[:, Hrep])
#     return (Hunique, Hrank)
# end

"""
    unique_haplotype_idx(H)

Returns the columns of `H` that are unique. 

# Input
* `H`: BitMatrix of haplotypes within a genomic window.

# Output
* BitVector where 1 indicates unique columns of H.
"""
function unique_haplotype_idx(H::BitMatrix)
    p = size(H, 1) 

    # reinterpret each haplotype as an integer
    if p == 8 
        HR = reinterpret(UInt8, H.chunks) 
    elseif p == 16
        HR = reinterpret(UInt16, H.chunks)
    elseif p == 32
        HR = reinterpret(UInt32, H.chunks)
    elseif p == 64
        HR = reinterpret(UInt64, H.chunks)
    elseif p == 128
        HR = reinterpret(UInt128, H.chunks)
    end

    return unique_index(HR)
end

function unique_index(v::AbstractVector)
    seen = Set{eltype(v)}()
    lv   = length(v)
    unique_index = trues(lv)

    @inbounds for i in 1:lv
        if in(v[i], seen)
            unique_index[i] = false
        else
            push!(seen, v[i])
        end
    end

    return unique_index
end


function redundant_haplotypes(
    X::AbstractMatrix{Union{Missing, T}},
    H::AbstractMatrix{T};
    width::Int    = 128,
    verbose::Bool = true
    ) where T <: Real

    # problem dimensions
    snps, people = size(X)

    # number of windows
    windows = ceil(Int, snps / width)

    # get unique haplotype indices and maps for each window
    Hunique  = unique_haplotypes(H, width, 'T')
    num_uniq = length(Hunique.uniqueindex[1])

    # Matrix storing redundant haplotypes. Each column is a person. Rows are redundant haplotypes for each window 
    redund_haps = PeoplesRedundantHaplotypeSet(windows, people) 

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
    compute_redundant_haplotypes!(redund_haps, Hunique, happair, H, 1)

    #TODO: make this loop multithreaded 
    for w in 2:(windows-1)

        # sync Xwork and Hwork with original data
        cur_range = ((w - 1) * width + 1):(w * width)
        M = resize_and_sync!(Xwork, Hwork, Hunique.uniqueindex[w], cur_range, X, H, M, N)

        # Calculate optimal haplotype pair among unique haplotypes
        haploimpute!(Xwork, Hwork, M, N, happair, hapscore, Xfloat=Xwork_float)

        # find all haplotypes matching the optimal haplotype pairs
        compute_redundant_haplotypes!(redund_haps, Hunique, happair, H, w)
    end

    # last window
    # TODO: there is problem with calling resize_and_sync! on the last window since the range is different
    last_range = ((windows - 1) * width + 1):snps
    M = resize_and_sync!(Xwork, Hwork, Hunique.uniqueindex[end], last_range, X, H, M, N)
    haploimpute!(Xwork, Hwork, M, N, happair, hapscore)
    compute_redundant_haplotypes!(redund_haps, Hunique, happair, H, windows)

    return redund_haps
end

# computational routine for recording redundant haplotypes for each window
function compute_redundant_haplotypes!(
    redund_haps::PeoplesRedundantHaplotypeSet, 
    Hunique::UniqueHaplotypeMaps, 
    happair::Tuple{AbstractVector, AbstractVector}, 
    H::AbstractMatrix,
    window::Int,
    )

    people = size(redund_haps, 2)

    # loop through all people
    @inbounds for k in 1:people
        (Hwork_i, Hwork_j) = (happair[1][k], happair[2][k])
        # println("person $k's optimal haplotype pairs are: $((Hwork_i, Hwork_j))")

        (H_i, H_j) = (Hunique.uniqueindex[window][Hwork_i], Hunique.uniqueindex[window][Hwork_j])
        # println("person $k's optimal haplotype pairs are located at columns $H_i and $H_j in H")

        # loop through all haplotypes and find ones that match either of the optimal haplotypes 
        for jj in 1:size(H, 2)
            Hunique.hapmap[window][jj] == H_i && push!(redund_haps.strand1[window, k], jj)
            Hunique.hapmap[window][jj] == H_j && push!(redund_haps.strand2[window, k], jj)
        end

        # println("person $k's redundant haplotypes are: ")
        # println(redund_haps[1, k])
    end

    return nothing
end

"""
    resize_and_sync!(X, H, M, N, Xwork, Hwork, Hnext, window)

Up/downsizes the dimension of `Hwork`, `M`, and `N` and copies relevant information into `Xwork` and `Hwork`. 

# Inputs
* `Xwork`: Worker matrix storing X[window, :]. 
* `Hwork`: Haplotype matrix in the current window containing only unique haplotypes. Must add/subtract columns. 
* `Hnext`: The unique haplotype indices of the current haplotype window. 
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
