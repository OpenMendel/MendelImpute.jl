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

function compute_optimal_halotype_set(
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
# function compute_redundant_haplotypes!(
#     redund_haps::PeoplesRedundantHaplotypeSet, 
#     Hunique::UniqueHaplotypeMaps, 
#     happair::Tuple{AbstractVector, AbstractVector}, 
#     H::AbstractMatrix,
#     window::Int,
#     )

#     people = size(redund_haps, 2)

#     # loop through all people
#     @inbounds for k in 1:people
#         (Hwork_i, Hwork_j) = (happair[1][k], happair[2][k])
#         # println("person $k's optimal haplotype pairs are: $((Hwork_i, Hwork_j))")

#         (H_i, H_j) = (Hunique.uniqueindex[window][Hwork_i], Hunique.uniqueindex[window][Hwork_j])
#         # println("person $k's optimal haplotype pairs are located at columns $H_i and $H_j in H")

#         # loop through all haplotypes and find ones that match either of the optimal haplotypes 
#         for jj in 1:size(H, 2)
#             Hunique.hapmap[window][jj] == H_i && push!(redund_haps.strand1[window, k], jj)
#             Hunique.hapmap[window][jj] == H_j && push!(redund_haps.strand2[window, k], jj)
#         end

#         # println("person $k's redundant haplotypes are: ")
#         # println(redund_haps[1, k])
#     end

#     return nothing
# end

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

# """
#     unique_haplotypes(H, window::UnitRange{Int})

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
