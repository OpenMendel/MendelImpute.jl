#testing hua's code

using Revise
using MendelImpute
using BenchmarkTools

n, p = 1000, 10000
X1 = rand(Float32, n, p);
X2 = Matrix{Union{Missing, eltype(X1)}}(X1);
X3 = ifelse.(rand(eltype(X1), n, p) .< 0.9, X2, missing);

@benchmark initmissing!(X3) setup=(X3=ifelse.(rand(n, p) .< 0.9, X2, missing))
@benchmark X4 = convert(Matrix{eltype(X1)}, X2)

using Revise
using BenchmarkTools
using MendelImpute
using Random
using LinearAlgebra

Random.seed!(123)
n = 5000 # number of individuals
p = 500  # number of SNPs
d = 500  # number of reference haplotypes
H = convert(Matrix{Float32}, rand(0:1, p, d))
X = convert(Matrix{Float32}, rand(0:2, p, n))
M = Transpose(H) * H
for j in 1:d, i in 1:(j - 1) # off-diagonal
    M[i, j] = 2M[i, j] + M[i, i] + M[j, j]
end
for j in 1:d # diagonal
    M[j, j] *= 4
end
N = Transpose(X) * H
for I in eachindex(N)
    N[I] *= 2
end
happair  = zeros(Int, n), zeros(Int, n)
hapscore = zeros(eltype(N), n)
@time haplopair!(happair, hapscore, M, N)

haplopair(X, H)




############ test best way to convert bitarray to BigInt

using BenchmarkTools
using Random

Random.seed!(123)
v = rand(0:1, 10)
u = rand(0:1, 512)

#add using BigInt
function bits_to_int(v::AbstractVector)
    l = length(v)
    total = BigInt(0)
    for i in 1:l
        cur = l - i + 1
        if v[i] == 1
            total += 2^BigInt(cur - 1)
        end
    end
    return total
end

#concatenate bits to string then parse to BigInt
function bits_to_int2(v::AbstractVector)
    concat = join(v)
    return parse(BigInt, concat, base=2)
end

#concatenate bits to string in buffer, then parse to BigInt
function bits_to_int3(v::AbstractVector)
    io = IOBuffer()
    for entry in v
        print(io, entry)
    end
    str = String(take!(io))
    concat = join(str)
    return parse(BigInt, concat, base=2)
end


bits_to_int(v) == bits_to_int2(v) == bits_to_int3(v) == 53
@benchmark bits_to_int(u) #5.305 ms, 9.78 MiB
@benchmark bits_to_int2(u) #1.074 ms, 957.39 KiB
@benchmark bits_to_int3(u) #1.244 ms, 975.97 KiB







#test best way to convert array of ints to bitarrays
using BenchmarkTools
using Random

Random.seed!(123)
z  = rand(1:10000, 10)
zz = rand(1:10000000, 1000000)

function int2bits(z::AbstractVector)
    uints = unique(z)
    nints = length(uints)
    uH    = zeros(eltype(z), 64, nints)
    for i in 1:nints
        uH[:, i] .= eltype(z).(digits(z[i], base=2, pad=64))
    end
    return uH
end

H = int2bits(z)



@benchmark int2bits(zz)



###### Fast Elimination of Redundant Haplotypes
using Revise
using BenchmarkTools
using MendelImpute
using Random
using LinearAlgebra
using Profile

#check correctness with eyeball
Random.seed!(13)
p = 8   # number of SNPs within a window
d = 5 # number of reference haplotypes
H = convert(Matrix{Float32}, rand(0:1, p, d))
H[:, 1] .= H[:, 2]
uH  = filter_redundant_haplotypes(H)
uH2 = filter_redundant_haplotypes(H)

#see if we can beat `unique` provided natively through julia
p = 512   # number of SNPs within a window
d = 10000 # number of reference haplotypes
H = convert(Matrix{Float32}, rand(0:1, p, d))
@benchmark unique(H, dims=1)               #65.839 ms , 19.57  MiB
@benchmark filter_redundant_haplotypes(H)  #631.095 ms, 533.01 MiB

p = 4   # number of SNPs within a window
d = 10000 # number of reference haplotypes
H = convert(Matrix{Float32}, rand(0:1, p, d))
@benchmark unique(H, dims=1)               #584.826 μs , 157.44 KiB
@benchmark filter_redundant_haplotypes(H)  #9.530 ms   , 8.46 MiB

p = 512   # number of SNPs within a window
d = 1_000_000 # number of reference haplotypes
H = convert(Matrix{Float32}, rand(0:1, p, d))
@benchmark unique(H, dims=1)               #7.756 s, 1.91 GiB
@benchmark filter_redundant_haplotypes(H)  #didn't finish in a few minutes






u = rand(1:5, 10)
v = rand(1:100, 10_000);
redundant_index(u)

function redundant_index(v::AbstractVector)
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

@benchmark redundant_index(v)




using Revise
using MendelImpute
using Random
using LinearAlgebra
using Profile
using StatsBase
using Random
using BenchmarkTools

Random.seed!(123)

function julia_unique(H)
    uH = unique(H, dims=1)
    return convert(Matrix{Float32}, uH)
end

p = 8     # number of SNPs within a window
d = 10000 # number of reference haplotypes
H = bitrand(p, d)
@benchmark unique_haplotypes(H)   #2.603 ms, 194.97 KiB
@benchmark julia_unique(H)        #2.196 ms, 323.75 KiB
@benchmark unique_haplotype_idx(H)#155.454 μs, 5.38 KiB

p = 64   # number of SNPs within a window
d = 10000 # number of reference haplotypes
H = bitrand(p, d)
@benchmark unique_haplotypes(H)   #9.545 ms, 3.77 MiB
@benchmark julia_unique(H)        #16.125 ms, 2.53 MiB
@benchmark unique_haplotype_idx(H)#445.094 μs, 195.28 KiB

p = 128   # number of SNPs within a window
d = 10000 # number of reference haplotypes
H = bitrand(p, d)
@benchmark unique_haplotypes(H)   #16.982 ms, 6.70 MiB
@benchmark julia_unique(H)        #34.475 ms, 5.05 MiB
@benchmark unique_haplotype_idx(H)#1.281 ms, 366.00 KiB







### search min in Matrix

using Revise
using BenchmarkTools
using MendelImpute
using Random
using LinearAlgebra
using Profile

Random.seed!(123)
n = 5000 # number of individuals
p = 500  # number of SNPs
d = 500  # number of reference haplotypes
H = convert(Matrix{Float32}, rand(0:1, p, d))
X = convert(Matrix{Float32}, rand(0:2, p, n))
M = Transpose(H) * H
for j in 1:d, i in 1:(j - 1) # off-diagonal
    M[i, j] = 2M[i, j] + M[i, i] + M[j, j]
end
for j in 1:d # diagonal
    M[j, j] *= 4
end
N = Transpose(X) * H
for I in eachindex(N)
    N[I] *= 2
end

happair  = zeros(Int, n), zeros(Int, n)
hapscore = zeros(eltype(N), n)
haplopair!(happair, hapscore, M, N)


M_col_min = zeros(eltype(N), d)
M_min_pos = zeros(Int, d)
happair2  = zeros(Int, n), zeros(Int, n)
hapscore2 = zeros(eltype(N), n)
haplopair2!(happair2, hapscore2, M, N, M_col_min=M_col_min, M_min_pos=M_min_pos)

all(happair .== happair2)
all(hapscore .== hapscore2)
[happair[1] happair[2] happair2[1] happair2[2]]
[hapscore hapscore2]

@benchmark haplopair!(happair, hapscore, M, N) #432.302 ms, 64 bytes
@benchmark haplopair2!(happair2, hapscore2, M, N, M_col_min=M_col_min, M_min_pos=M_min_pos) #542.754 ms, 64 bytes





function run()
    b = ones(100)
    for i in 1:100
        x .= U \ (L \ b)
    end
end





using Revise
using BenchmarkTools
using MendelImpute
using Random
using LinearAlgebra
using Profile

Random.seed!(123)
n, p, d = 5000, 500, 500
H = convert(Matrix{Float32}, rand(0:1, p, d))
X = convert(Matrix{Float32}, rand(0:2, p, n))
M = zeros(eltype(H), d, d)
N = zeros(promote_type(eltype(H), eltype(X)), n, d)
happair  = zeros(Int, n), zeros(Int, n)
hapscore = zeros(eltype(N), n)
missingprop = 0.1

X2 = Matrix{Union{Missing, eltype(X)}}(X)
X3 = ifelse.(rand(eltype(X), p, n) .< missingprop, missing, X2)
X3_original = deepcopy(X3)

missing_location = BitArray(undef, size(X))
findmissing!(X3, missing_location)

@benchmark findmissing!(X3, missing_location) setup=(missing_location = BitArray(undef, size(X3)))






using Revise
using BenchmarkTools
using MendelImpute
using Random
using LinearAlgebra
using Profile

Random.seed!(123)
d = 10000
H = bitrand(128000, d)
# Hwork = unique_haplotypes(H, 1:128)
# Hunique = unique_haplotypes(H, 128)

for i in 1:Int(d/2)
    H[:, 2i] .= H[:, 2i - 1]
end

windows = 10
width = 64
storage = BitMatrix(undef, 64, size(H, 2))
copyto!(storage, @view(H[1:width, :]))


p, d    = size(H)
windows = ceil(Int, p / width)
unique_hap = UniqueHaplotypes(windows, d)
fast_data_type = Dict(8=>UInt8, 16=>UInt16, 32=>UInt32, 64=>UInt64, 128=>UInt128)

HR = reinterpret(fast_data_type[width], storage.chunks)


#128000 by 10000 H, width = 64,  43.702105 seconds (47.99 k allocations: 378.702 MiB, 0.57% gc time)
#128000 by 10000 H, width = 128, 35.878729 seconds (49.97 k allocations: 712.912 MiB, 0.96% gc time)
@time unique_hap = fast_elimination(H, windows, width, H[1:width, :], fast_data_type)

#128000 by 10000 H, width = 64, 15.323239 seconds (272.48 k allocations: 3.146 GiB, 3.46% gc time)
#128000 by 10000 H, width = 128, 14.234709 seconds (135.48 k allocations: 1.573 GiB, 2.94% gc time)
@time Hunique = unique_haplotypes(H, 64, 'T')

# Hwork = unique_haplotypes(H, 1:128, 'T')
Hunique = unique_haplotypes(H, 64, 'T')

storage = zeros(Int, 1000)
hii = groupslices(H, 2)
groupslices!(storage, H, 2)
# unique(groupslices(H, 2))
# unique_haplotype_idx(H)

@benchmark unique_haplotypes(H, 1:128)

@benchmark unique(groupslices(H, 2)) # 1.273 ms, 348.78 KiB
@benchmark unique(H, dims=2)         # 1.408 ms, 132.25 KiB
@benchmark unique_haplotype_idx(H)   # 88.825 μs, 92.58 KiB




H = zeros(1, 9)
H[1] = 1
H[2] = 2
H[3] = 2
H[4] = 3
H[5] = 2
H[6] = 4
H[7] = 1
H[8] = 4
H[9] = 1

result = unique_haplotypes(H, 128, 'T')




#AFRped data

using Revise
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using MendelImpute
using Random
using Profile
using ElasticArrays


rawdata = readdlm("AFRped_geno.txt", ',', Float32);
people = 664;
X = copy(Transpose(rawdata[1:people, 1:(end - 1)]));
function create_hap(x)
    n, p = size(x)
    h = one(eltype(x))
    for j in 1:p, i in 1:n
        if x[i, j] != 0
            x[i, j] -= h
        end
    end
    return copy(Transpose(x))
end
H = create_hap(rawdata[(people + 1):end, 1:(end - 1)]);

Random.seed!(123)
missingprop = 0.1
p, n = size(X)
X2 = Matrix{Union{Missing, eltype(X)}}(X)
Xm = ifelse.(rand(eltype(X), p, n) .< missingprop, missing, X2)
Xm_original = copy(Xm)

# Hunique = unique_haplotypes(H, 128, 'T')
# @benchmark unique_haplotypes(H, 128, 'T') #247.271 ms, 16.24 MiB
# @benchmark unique_haplotypes(H, 64, 'T')  #240.313 ms, 20.95 MiB
# @benchmark unique_haplotypes(H, 32, 'T')  #275.781 ms, 29.25 MiB
# @benchmark unique_haplotypes(H, 16, 'T')  #287.333 ms, 48.92 MiB

ph2 = phase2(Xm, H, width=32)

@time ph2 = phase2(Xm, H, width=32, verbose=false); # downsizing M is amortized: 1.096796 seconds (3.92 M allocations: 176.796 MiB)
@time ph2 = phase2(Xm, H, width=32, verbose=false); # always reallocate new M:   1.102165 seconds (3.92 M allocations: 177.990 MiB)
@time ph2 = phase2(Xm, H, width=32, verbose=false); # calls resize! on Mvec:     1.095878 seconds (3.92 M allocations: 174.128 MiB)

@time ph2 = phase2(Xm, H, width=128, verbose=false); # downsizing M is amortized:   1.403395 seconds (797.02 k allocations: 77.320 MiB)
@time ph2 = phase2(Xm, H, width=128, verbose=false); # always reallocate new M:     1.360186 seconds (796.94 k allocations: 85.620 MiB)
@time ph2 = phase2(Xm, H, width=128, verbose=false); # calls resize! on Mvec:       1.392999 seconds (796.71 k allocations: 57.788 MiB)


#resize matrix efficiently

using BenchmarkTools
using LinearAlgebra
using ElasticArrays

H = rand(1000, 1000)
Hwork = ElasticArray{Float64}(H[1:100, 1:100])
resize!(Hwork, 100, 1000)





A = rand(1000, 1000)
H = zeros(1000, 1000)
Htest = ElasticArray{Float64}(undef, 1000, 1000)
@benchmark copyto!(H, A)
@benchmark copyto!(Htest, A)

M = rand(100, 100)
Mvec = vec(M)
resize!(Mvec, 1000000)
M = Base.ReshapedArray(Mvec, (1000, 1000), ())

C = zeros(1000, 1000)
A = rand(1000, 1000)
B = rand(1000, 1000)
@benchmark mul!(C, A, B) #9.462 ms
@benchmark mul!(C, M, B) #10.058 ms





using Revise
using MendelImpute
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using Random
using ElasticArrays

cd("/Users/biona001/.julia/dev/MendelImpute/data")

rawdata = readdlm("AFRped_geno.txt", ',', Float32);
people = 664;
X = copy(Transpose(rawdata[1:people, 1:(end - 1)]));
function create_hap(x)
    n, p = size(x)
    h = one(eltype(x))
    for j in 1:p, i in 1:n
        if x[i, j] != 0
            x[i, j] -= h
        end
    end
    return copy(Transpose(x))
end
H = create_hap(rawdata[(people + 1):end, 1:(end - 1)]);

Random.seed!(123)
missingprop = 0.1
p, n = size(X)
X2 = Matrix{Union{Missing, eltype(X)}}(X)
Xm = ifelse.(rand(eltype(X), p, n) .< missingprop, missing, X2)
Xm_original = copy(Xm)
width = 64
windows = floor(Int, p / width)

# computes optimal-redundant haplotypes for each window/person
@time opt = compute_optimal_halotype_set(Xm, H, width=width, verbose=true) #2.187847 seconds (1.55 M allocations: 184.291 MiB)
opt[1].strand1 #person 1's optimal haplotypes on strand1 for each window

#Hua's code error = 0.013899062884015356 (width = 64)
#Hua's code error = 0.0033518189729195617 (width = 400) #12.170407 seconds (87.59 k allocations: 23.879 MiB, 0.15% gc time)
@time ph = phase(Xm, H, width)

hapset, phase = phase2(Xm, H, width=64);
hapset, phase = phase2(Xm, H, width=3*64);
hapset, phase = phase2(Xm, H, width=400);
hapset, phase = phase2(Xm, H, width=1200);

@benchmark phase2(Xm, H, width=64) seconds=15   # width 64  : 2.925 s, 215.08 MiB, 1941655 alloc
@benchmark phase2(Xm, H, width=400) seconds=15  # width 400 : 5.149 s, 62.97 MiB, 390094 alloc
@benchmark phase2(Xm, H, width=1200) seconds=15 # width 1200: 4.603 s, 53.94 MiB, 153625 alloc

# look at the haplotype intersections
findfirst.(hapset[1].strand1)
findfirst.(hapset[1].strand2)

findfirst.(hapset[10].strand1)
findfirst.(hapset[10].strand2)

ph = phase(Xm, H, width=width)

# calculate error rate
impute2!(Xm, H, ph)
missing_idx    = ismissing.(Xm_original)
total_missing  = sum(missing_idx)
actual_missing_values  = convert(Vector{Int64}, X[missing_idx])  #true values of missing entries
imputed_missing_values = convert(Vector{Int64}, Xm[missing_idx]) #imputed values of missing entries
error_rate = sum(actual_missing_values .!= imputed_missing_values) / total_missing
copyto!(Xm, Xm_original);

# profiling
phase2(Xm, H, width=width)
using Profile
Profile.clear()
@profile phase2(Xm, H, width=width)
using ProfileView
ProfileView.view()

# old code error:
#current code error = 0.019282749489147502 (width = 64)
#current code error = 0.008566492445731325 (width = 3*64)
#current code error = 0.0045481021680262605 (width = 400)
#current code error = 0.003721998955416397 (width = 3*400) #3.769867 seconds (663.78 k allocations: 61.506 MiB, 0.80% gc time)

# Without searching for breakpoints:
# error = 0.030189043249636106 (width = 64) ====> this should get down to 0.026 since that's the error in old code
# error = 0.014054060293167706 (width = 3*64)
# error = 0.007871065240305758 (width = 400)
# error = 0.0063462370050543174 (width = 1200)

# Searching for breakpoints:
# error = 0.019286459533515512 (width = 64) ====> this should get down to 0.0192 since that's the error in old code
# error = 0.008584218213267367 (width = 3*64)
# error = 0.004584790384554343 (width = 400)
# error = 0.0038312391506966433 (width = 1200)

function naive_impute!(X)
    n, p = size(X)
    fillval = convert(eltype(X), 0.0)
    for j in 1:p, i in 1:n
        ismissing(X[i, j]) && (X[i, j] = fillval)
    end
end
naive_impute!(Xm) #fill with 0 gives error rate = 0.09844272948788609
naive_impute!(Xm) #fill with 1 gives error rate = 0.9343025343313078
naive_impute!(Xm) #fill with 2 gives error rate = 0.9672547361808062



function naive_impute!(X, X_original)
    n, p = size(X)
    fillval = 0x00
    for j in 1:p, i in 1:n
        if X_original[i, j] == 0x01
            X[i, j] = fillval
        else
            X[i, j] = X_original[i, j]
        end
    end
end




a = BitVector(undef, 1_000_000)
b = bitrand(1_000_000)
c = bitrand(1_000_000)
@benchmark $a .= $b .& $c

function test(H::BitMatrix)
    println("H is BitMatrix")
end

function test(H::AbstractMatrix)
    println("H is AbstractMatrix")
end



# take away: don't use sparse arrays
using SparseArrays
using BenchmarkTools

n = 50_000
a = falses(n)
b = falses(n)
c = falses(n)

b[1:100] .= true
c[1:100] .= true
shuffle!(b)
shuffle!(c)

sa = sparse(a)
sb = sparse(b)
sc = sparse(c)

@benchmark $sa .= $sb .& $sc
@benchmark $a .= $b .& $c


using BenchmarkTools

n = 50_000
a = BitVector(undef, n);
b = rand(1:n, n);
c = 1

@benchmark $a .= $b .== $c




# simulate utilities test
using Revise
using MendelImpute
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using Random
using ElasticArrays

cd("/Users/biona001/.julia/dev/MendelImpute/data")

rawdata = readdlm("AFRped_geno.txt", ',', Float32);
people = 665;
X = copy(Transpose(rawdata[1:people, 1:(end - 1)]));
function create_hap(x)
    n, p = size(x)
    h = one(eltype(x))
    for j in 1:p, i in 1:n
        if x[i, j] != 0
            x[i, j] -= h
        end
    end
    return copy(Transpose(x))
end
H = create_hap(rawdata[(people + 1):end, 1:(end - 1)]);

# tgt = convert(Matrix{Int}, X)
# ref = convert(BitArray{2}, H)
# make_refvcf_file(ref, filename="AFRped_ref.vcf")
# make_tgtvcf_file(tgt, filename="AFRped_tgt.vcf")

Random.seed!(123)
missingprop = 0.1
p, n = size(X)
X2 = Matrix{Union{Missing, eltype(X)}}(X)
Xm = ifelse.(rand(eltype(X), p, n) .< missingprop, missing, X2)
Xm_original = copy(Xm)
width = 64
windows = floor(Int, p / width)

hapset, phase = phase2(Xm, H, width=64);  #error = 0.019300418, 215.47MB, 3.070s, 1944069 alloc
hapset, phase = phase2(Xm, H, width=700); #error = 0.003373971, 48.53MB, 4.880s, 205942 alloc

# calculate error rate
impute2!(Xm, H, phase)
missing_idx    = ismissing.(Xm_original)
total_missing  = sum(missing_idx)
actual_missing_values  = convert(Vector{Int64}, X[missing_idx])  #true values of missing entries
imputed_missing_values = convert(Vector{Int64}, Xm[missing_idx]) #imputed values of missing entries
error_rate = sum(actual_missing_values .!= imputed_missing_values) / total_missing
copyto!(Xm, Xm_original);


#simulation code
using Revise
using MendelImpute
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using Random
using ElasticArrays

p = 10000
d = 250

Random.seed!(2019)

H = simulate_markov_haplotypes(p, d)
X = simulate_genotypes(H, 100)

count(iszero, H) #124388
count(isone, H) #125612

count(iszero, X) #24774
count(isone, X) #50058
count(x -> x == 2, X) #25168

Random.seed!(2019)
H = simulate_uniform_haplotypes(p, d)
X = simulate_genotypes(H, people = 100)

count(iszero, H) #187795
count(isone, H) #62205

count(iszero, X) #56689
count(isone, X) #37049
count(x -> x == 2, X) #6262




## test writer
using Revise
using MendelImpute
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using Random
using ElasticArrays
using GeneticVariation
using VCFTools

cd("/Users/biona001/.julia/dev/MendelImpute/data")

rawdata = readdlm("AFRped_geno.txt", ',', Float32);
people = 665;
X = copy(Transpose(rawdata[1:people, 1:(end - 1)]));
function create_hap(x)
    n, p = size(x)
    h = one(eltype(x))
    for j in 1:p, i in 1:n
        if x[i, j] != 0
            x[i, j] -= h
        end
    end
    return copy(Transpose(x))
end
H = create_hap(rawdata[(people + 1):end, 1:(end - 1)]);

Random.seed!(123)
missingprop = 0.1
p, n = size(X)
X2 = Matrix{Union{Missing, eltype(X)}}(X)
Xm = ifelse.(rand(eltype(X), p, n) .< missingprop, missing, X2)
Xm_original = copy(Xm)
width = 64
windows = floor(Int, p / width)

hapset, phase = phase2(Xm, H, width=700);

vcffile = "test.08Jun17.d8b.vcf"
des = "phased." * vcffile
reader = VCF.Reader(openvcf(vcffile, "r"))
writer = VCF.Writer(openvcf(des, "w"), reader.header)

@time hapset, phase = phase2(Xm, H, width=700); #4.497113 seconds (205.95 k allocations: 48.530 MiB, 0.38% gc time)
@time write(writer, phase, H) #7.760277 seconds (97.08 M allocations: 4.340 GiB, 5.29% gc time)



io = open("myfile.txt", "w");
x = 0
write(io, x, '\n')
print(io, x, '\t', 0.0, "hi");
close(io);




using Revise
using MendelImpute
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using Random
using ElasticArrays

cd("/Users/biona001/.julia/dev/MendelImpute/data")

rawdata = readdlm("AFRped_geno.txt", ',', Float32);
people = 664;
X = copy(Transpose(rawdata[1:people, 1:(end - 1)]));
function create_hap(x)
    n, p = size(x)
    h = one(eltype(x))
    for j in 1:p, i in 1:n
        if x[i, j] != 0
            x[i, j] -= h
        end
    end
    return copy(Transpose(x))
end
H = create_hap(rawdata[(people + 1):end, 1:(end - 1)]);

Random.seed!(123)
missingprop = 0.1
p, n = size(X)
X2 = Matrix{Union{Missing, eltype(X)}}(X)
Xm = ifelse.(rand(eltype(X), p, n) .< missingprop, missing, X2)
Xm_original = copy(Xm)
width = 64
windows = floor(Int, p / width)

ph = phase(Xm, H, width=64)

@benchmark phase(Xm, H, width=64) seconds=15   # width 64  : 3.914 s, 401.31 MiB, 11887397 alloc (this includes calculating optimal hapset)
@benchmark phase(Xm, H, width=400) seconds=15  # width 400 :
@benchmark phase(Xm, H, width=1200) seconds=15 # width 1200:


using Revise
using MendelImpute
H = simulate_markov_haplotypes(10000, 500)
X = simulate_genotypes(H, 100)
X2 = simulate_genotypes2(H, people=100)



using Revise
using MendelImpute
using VCFTools
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using Random
using ElasticArrays

cd("/Users/biona001/.julia/dev/MendelImpute/data")

file = "test.08Jun17.d8b.vcf.gz"
X = convert_ht(Float32, file, trans=true)
H = convert_ht(Float32, file, trans=true)

width = 400
hapset = compute_optimal_halotype_set(X, H, width=width)



# simulate haplotypes in windows for visualization

using Revise
using MendelImpute
using VCFTools
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using Random
using ElasticArrays
using StatsBase

Random.seed!(2020)
windows = 10

# generate happairs in windows
T = Tuple{Int, Int}
windows = 5
haplotype_set = Vector{Vector{T}}(undef, windows)
pu

# strand 1
s1 = [Int[] for w in 1:windows]
for w in 1:windows
    nhaps = rand(2:5)
    haps  = sample(1:10, nhaps, replace=false) # sample 2~5 haplotypes from 1:10
    sort!(haps)
    for i in haps
        push!(s1[w], i)
    end
end
s1

# strand 2
s2 = [Int[] for w in 1:windows]
for w in 1:windows
    nhaps = rand(2:5)
    haps  = sample(1:10, nhaps, replace=false) # sample 2~5 haplotypes from 1:10
    sort!(haps)
    for i in haps
        push!(s2[w], i)
    end
end
s2




# simulate haplotypes in windows for visualization

using Revise
using MendelImpute
using VCFTools
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using Random
using ElasticArrays
using StatsBase

# generate happairs in windows
T = Tuple{Int, Int}
windows = 4
haplotype_set = [T[] for i in 1:windows]
push!(haplotype_set[1], (1, 4))
push!(haplotype_set[1], (1, 5))
push!(haplotype_set[1], (1, 6))
push!(haplotype_set[1], (2, 4))
push!(haplotype_set[1], (2, 5))
push!(haplotype_set[1], (2, 6))
push!(haplotype_set[1], (3, 4))
push!(haplotype_set[1], (3, 5))
push!(haplotype_set[1], (3, 6))

push!(haplotype_set[2], (1, 5))
push!(haplotype_set[2], (1, 7))
push!(haplotype_set[2], (1, 8))
push!(haplotype_set[2], (2, 5))
push!(haplotype_set[2], (2, 7))
push!(haplotype_set[2], (2, 8))
push!(haplotype_set[2], (6, 5))
push!(haplotype_set[2], (6, 7))
push!(haplotype_set[2], (6, 8))

push!(haplotype_set[3], (1, 1))
push!(haplotype_set[3], (1, 5))
push!(haplotype_set[3], (3, 1))
push!(haplotype_set[3], (3, 5))

push!(haplotype_set[4], (4, 4))
push!(haplotype_set[4], (4, 6))
push!(haplotype_set[4], (4, 8))
push!(haplotype_set[4], (5, 4))
push!(haplotype_set[4], (5, 6))
push!(haplotype_set[4], (5, 8))
push!(haplotype_set[4], (8, 4))
push!(haplotype_set[4], (8, 6))
push!(haplotype_set[4], (8, 8))
haplotype_set

sol_path, next_pair, subtree_err, best_err = connect_happairs(haplotype_set)



using Revise
using MendelImpute
using VCFTools
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using Random
using ElasticArrays
using StatsBase

# generate happairs in windows
T = Tuple{Int, Int}
windows = 5
haplotype_set = [T[] for i in 1:windows]

Random.seed!(2020)
for w in 1:windows
    haplotype_set[w] = [(rand(1:10), rand(1:10)) for i in 1:rand(1:10)]
end
sol_path, next_pair, subtree_err, best_err = connect_happairs(haplotype_set)



# generate happairs in windows
windows = 10
T = Tuple{Int, Int}
haplotype_set = [T[] for i in 1:windows]

Random.seed!(123)
for w in 1:windows
    haplotype_set[w] = [(rand(1:10), rand(1:10)) for i in 1:rand(1:10)]
end
haplotype_set
sol_path, next_pair, subtree_err, best_err = connect_happairs(haplotype_set)

# try another 10window test
windows = 10
T = Tuple{Int, Int}
haplotype_set = [T[] for i in 1:windows]
Random.seed!(321)
for w in 1:windows
    haplotype_set[w] = [(rand(1:10), rand(1:10)) for i in 1:rand(1:10)]
end
haplotype_set
sol_path, memory, best_err = connect_happairs(haplotype_set)



# generate happairs in windows
windows = 10000
haplotype_set = [T[] for i in 1:windows]

Random.seed!(2020)
for w in 1:windows
    haplotype_set[w] = [(rand(1:10000), rand(1:10000)) for i in 1:rand(1:100)]
end
haplotype_set

sol_path1, memory, path_err, best_err = connect_happairs(haplotype_set)#636.229 ms, 61.76 MiB
sol_path2, memory, best_err = connect_happairs2(haplotype_set) #333.806 ms, 95.85 MiB
sum(sol_path1 .!= sol_path2) # 89 places different

# vector of sets vs vector of vector of tuples
using BenchmarkTools
n = 1000
T = Tuple{Int, Int}
@benchmark x = [Set{T}() for i in 1:n] # 4.291 ms, 617.44 KiB
@benchmark y = [T[] for i in 1:n]      # 72.620 μs, 86.19 KiB
@benchmark push!($x[1], (1, 2)) # 11.392 ns, 0 bytes
@benchmark push!($y[1], (1, 2)) # 15.803 ns, 0 bytes


function test_push(x)
    for i in 1:1000
        push!(x[rand(1:1000)], (rand(1:100000), rand(1:100000)))
    end
    return x
end
x = [Set{T}() for i in 1:n]
y = [T[] for i in 1:n]
@benchmark test_push(x) # 205.370 μs
@benchmark test_push(y) # 47.175 μs




using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
cd("/Users/biona001/.julia/dev/MendelImpute/simulation")

# impute
tgtfile = "./compare6/target_masked.vcf.gz"
reffile = "./compare6/haplo_ref.vcf"
outfile = "./compare6/imputed_target.vcf.gz"
width   = 800

H = convert_ht(Float64, reffile, trans=true)
X = convert_gt(Float64, tgtfile, trans=true)
hapset = compute_optimal_halotype_set(X, H, width = width)
happair = compute_optimal_halotype_pair(X, H, width = width)

# hapset[1][1]
@time sol_path, memory, path_err, best_err = connect_happairs(hapset[1]); # 0.077636 seconds (201 allocations: 676.344 KiB)
@time sol_path, memory, best_err = connect_happairs2(hapset[1]); # 0.039125 seconds (209 allocations: 1.071 MiB)


@time hs, ph = phase(tgtfile, reffile, impute=true, outfile = outfile, width = width);

# import imputed result and compare with true
X_complete  = convert_gt(Float32, "./compare6/target.vcf.gz");
X_mendel = convert_gt(Float32, outfile);
n, p = size(X_mendel)
error_rate = sum(X_mendel .!= X_complete) / n / p

# not searching breakpoints (w = 400):
# 3462.416399 seconds (77.81 G allocations: 3.402 TiB, 12.66% gc time)
# error = 0.00037152087475149104

# searching only 1 strand's bkpt (w = 800):
# 35.990546 seconds (98.16 M allocations: 9.233 GiB, 3.80% gc time)
# error = 0.0005958889520022721

# searching both strand's bkpt (w = 800)
# 167.675135 seconds (98.00 M allocations: 9.224 GiB, 0.75% gc time)
# error = 0.0002344859414938938



# different haplopair! strategy:
# keep best pair only (orignal code): error = 0.00025205907412666857, 187.475166 sec
# keep all happairs that are equally good: error = 0.0002509940357852883, 183.600139 sec
# keep top 10 haplotype pairs: 0.0003010508378301619, 175.682698 sec
# keep all previous best pairs: error = 0.00023839108207895484, 186.741630 sec
# keep all previous best pairs and equally good pairs: 0.0002378585629082647, 189.419958 sec


using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
cd("/Users/biona001/.julia/dev/MendelImpute/simulation")

# impute
tgtfile = "./compare6/target_masked.vcf.gz"
reffile = "./compare6/haplo_ref.vcf"
outfile = "./compare6/imputed_target.vcf.gz"
width   = 800

@time hs, ph = phase(tgtfile, reffile, impute=true, outfile = outfile, width = width);

# import imputed result and compare with true
X_complete  = convert_gt(Float32, "./compare6/target.vcf.gz"; as_minorallele=false);
X_mendel = convert_gt(Float32, outfile, as_minorallele=false);
n, p = size(X_mendel);
error_rate = sum(X_mendel .!= X_complete) / n / p

# print one person's error in each window
windows = floor(Int, p / width)
person = 3
for w in 1:windows
    win_range = ((w - 1) * width + 1):(w * width)
    error_rate = sum(X_complete[person, win_range] .!= X_mendel[person, win_range]) / length(win_range)
    println("person $person window $w = $win_range has error = $error_rate")
end
tot_error = sum(X_complete[person, :] .!= X_mendel[person, :]) / p

# calculate error for 1 window in 1 person skipmissing
x = X[5601:6400, 1]
h1 = H[5601:6400, 1165]
h2 = H[5601:6400, 542]
function run()
    errors = 0
    for pos in 1:n
        if !ismissing(x[pos])
            errors += x[pos] ≠ h1[pos] + h2[pos]
        end
    end
    println(errors)
end
run()

# print error for everybody
for person in 1:people
    error_rate = sum(X_complete[person, :] .!= X_mendel[person, :]) / p
    println("person $person error = $error_rate")
end




######## try to make dynamic programming efficient

# profile mamory usage
julia --track-allocation=user

using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using Profile
using BenchmarkTools

cd("/Users/biona001/.julia/dev/MendelImpute/simulation")
tgtfile = "./compare6/target_masked.vcf.gz"
reffile = "./compare6/haplo_ref.vcf"
outfile = "./compare6/imputed_target.vcf.gz"
width   = 800

H = convert_ht(Float64, reffile, trans=true)
X = convert_gt(Float64, tgtfile, trans=true)
hapset = compute_optimal_halotype_set(X, H, width = width)

T = Tuple{Int, Int}
windows  = length(hapset[3])
sol_path = Vector{Tuple{Int, Int}}(undef, windows)
next_pair = [Int[] for i in 1:windows]
subtree_err = [Float64[] for i in 1:windows]
@time MendelImpute.connect_happairs!(sol_path, next_pair, subtree_err, hapset[3]);
# 0.019480 seconds (7 allocations: 624 bytes)

@benchmark MendelImpute.connect_happairs!(sol_path, next_pair, subtree_err, hapset[3])
#   memory estimate:  464 bytes
#   allocs estimate:  3
#   --------------
#   minimum time:     11.387 ms (0.00% GC)
#   median time:      12.320 ms (0.00% GC)
#   mean time:        12.805 ms (0.00% GC)
#   maximum time:     28.351 ms (0.00% GC)
#   --------------
#   samples:          391
#   evals/sample:     1

Profile.clear_malloc_data()
MendelImpute.connect_happairs!(hapset[3], memory=mymemory, sol_path=sol_path, path_err=path_err);



# test if adding/looking up dictionaries is efficient
function mytest()
    mydict = [Dict{Tuple{Int, Int}, Float64}() for i in 1:length(hapset[3])]

    my_arbitrary_sum = 0.0
    for i in 1:length(hapset[3])
        # add to dictionaries
        num = length(hapset[3][i])
        for n in 1:num
            mydict[i][(rand(1:100), rand(1:100))] = rand() # add arbitrary pair
        end

        # lookup a bunch of times
        for n in 1:num
            pair = (rand(1:100), rand(1:100))
            if haskey(mydict[i], pair)
                my_arbitrary_sum += pair[1] + pair[2]
            end
        end
    end

    return mydict, my_arbitrary_sum
end
mydict, my_arbitrary_sum = mytest()
@benchmark mytest() # 2.728 ms, 2.58 MiB



# profile mamory usage
julia --track-allocation=user

using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using BenchmarkTools
using Random
using Profile

cd("/Users/biona001/.julia/dev/MendelImpute/simulation")
tgtfile = "./compare3/target_masked.vcf.gz"
reffile = "./compare3/haplo_ref.vcf.gz"
outfile = "./compare3/imputed_target.vcf.gz"
width   = 400

@time hs, ph = phase(tgtfile, reffile, impute=true, outfile = outfile, width = width);

Profile.clear_malloc_data()
hs, ph = phase(tgtfile, reffile, impute=true, outfile = outfile, width = width);





using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using Profile

cd("/Users/biona001/.julia/dev/MendelImpute/simulation")
tgtfile = "./compare6/target_masked.vcf.gz"
reffile = "./compare6/haplo_ref.vcf"
outfile = "./compare6/imputed_target.vcf.gz"
width   = 800

H = convert_ht(Float64, reffile, trans=true)
X = convert_gt(Float64, tgtfile, trans=true)
hapset = compute_optimal_halotype_set(X, H, width = width);

@time ph = phase(X, H, hapset=hapset, width=width);
# 128.425430 seconds (1.34 M allocations: 69.787 MiB)

@time hs, ph = phase(tgtfile, reffile, impute=true, outfile = outfile, width = width);




# checking memory allocations
using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using Profile

cd("/Users/biona001/.julia/dev/MendelImpute/simulation")
tgtfile = "./compare3/target_masked.vcf.gz"
reffile = "./compare3/haplo_ref.vcf.gz"
outfile = "./compare3/imputed_target.vcf.gz"
width   = 800

@time hs, ph = phase(tgtfile, reffile, impute=true, outfile = outfile, width = width);
# overall: 4.610923 seconds (5.37 M allocations: 517.820 MiB, 1.34% gc time)

@time H = convert_ht(Float32, reffile, trans=true); # 0.271169 seconds (3.07 M allocations: 285.923 MiB, 15.66% gc time)
@time X = convert_gt(Float32, tgtfile, trans=true); # 0.081043 seconds (1.26 M allocations: 114.565 MiB)
@time hapset = compute_optimal_halotype_set(X, H, width = width); # 0.318205 seconds (16.74 k allocations: 56.280 MiB)
@time ph = phase(X, H, hapset = hapset, width = width); # 4.097682 seconds (1.93 k allocations: 118.406 KiB)
@time write_test(ph, H) # 0.256324 seconds (941.18 k allocations: 99.117 MiB, 7.53% gc time)

function write_test(ph, H)
    reader = VCF.Reader(openvcf(tgtfile, "r"))
    writer = VCF.Writer(openvcf(outfile, "w"), header(reader))

    # loop over each record
    for (i, record) in enumerate(reader)
        gtkey = VCF.findgenokey(record, "GT")
        if !isnothing(gtkey)
            # loop over samples
            for (j, geno) in enumerate(record.genotype)
                # if missing = '.' = 0x2e
                if record.data[geno[gtkey][1]] == 0x2e
                    #find where snp is located in phase
                    hap1_position = searchsortedlast(ph[j].strand1.start, i)
                    hap2_position = searchsortedlast(ph[j].strand2.start, i)

                    #find the correct haplotypes
                    hap1 = ph[j].strand1.haplotypelabel[hap1_position]
                    hap2 = ph[j].strand2.haplotypelabel[hap2_position]

                    # save actual allele to data. "0" (REF) => 0x30, "1" (ALT) => 0x31
                    a1, a2 = convert(Bool, H[i, hap1]), convert(Bool, H[i, hap2])
                    record.data[geno[gtkey][1]] = ifelse(a1, 0x31, 0x30)
                    record.data[geno[gtkey][2]] = 0x7c # phased data has separator '|'
                    record.data[geno[gtkey][3]] = ifelse(a2, 0x31, 0x30)
                end
            end
        end
        write(writer, record)
    end

    # close
    flush(writer); close(reader); close(writer)
end



julia --track-allocation=user

using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using Profile

cd("/Users/biona001/.julia/dev/MendelImpute/simulation")
tgtfile = "./compare3/target_masked.vcf.gz"
reffile = "./compare3/haplo_ref.vcf.gz"
outfile = "./compare3/imputed_target.vcf.gz"
width   = 800


@time H = convert_ht(Float64, reffile, trans=true); # 0.271169 seconds (3.07 M allocations: 285.923 MiB, 15.66% gc time)
@time X = convert_gt(Float64, tgtfile, trans=true); # 0.081043 seconds (1.26 M allocations: 114.565 MiB)

@time hapset = compute_optimal_halotype_set(X, H, width = width); # 0.318205 seconds (16.74 k allocations: 56.280 MiB)
@time ph = phase(X, H, hapset = hapset, width = width); # 3.662491 seconds (1.93 k allocations: 118.406 KiB)

Profile.clear_malloc_data()
ph = phase(X, H, hapset = hapset, width = width); # 3.662491 seconds (1.93 k allocations: 118.406 KiB)




# optimize search_breakpoint(X, H, s1::Tuple{Int, Int}, s2::Tuple{Int, Int})
julia --track-allocation=user

using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using Profile
using BenchmarkTools

cd("/Users/biona001/.julia/dev/MendelImpute/simulation")
tgtfile = "./compare3/target_masked.vcf.gz"
reffile = "./compare3/haplo_ref.vcf.gz"
outfile = "./compare3/imputed_target.vcf.gz"
width   = 800

@time H = convert_ht(Float64, reffile, trans=true); # 0.271169 seconds (3.07 M allocations: 285.923 MiB, 15.66% gc time)
@time X = convert_gt(Float64, tgtfile, trans=true); # 0.081043 seconds (1.26 M allocations: 114.565 MiB)

Xi = @view(X[1:2width, 1])
Hi = @view(H[1:2width, :])
@code_warntype search_breakpoint(Xi, Hi, (1, 3), (4, 5))
@code_warntype search_breakpoint(Xi, Hi, 1, (4, 5))

@time search_breakpoint(Xi, Hi, (1, 3), (4, 5)) # 0.018376 seconds (7 allocations: 256 bytes)
@time search_breakpoint(Xi, Hi, 1, (4, 5))      # 0.000036 seconds (6 allocations: 224 bytes)

@benchmark search_breakpoint(Xi, Hi, (1, 3), (4, 5)) # 11.500 ms, 32 bytes, 1 alloc


# measure allocation
search_breakpoint(Xi, Hi, (1, 3), (4, 5))
Profile.clear_malloc_data()
search_breakpoint(Xi, Hi, (1, 3), (4, 5))



function change_tuples()
    x = (1, 2)
    for i in 1:100
        y = (rand(1:10000), rand(1:10000))
        x = y
    end
    return x
end

@time change_tuples()






using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using Profile

cd("/Users/biona001/.julia/dev/MendelImpute/simulation")
tgtfile = "target_masked.vcf.gz"
reffile = "haplo_ref.vcf.gz"
outfile = "imputed_target.vcf.gz"
width   = 400


# overall
@time hs, ph = phase(tgtfile, reffile, outfile = outfile, width = width);
# 112.814546 seconds (368.18 M allocations: 34.114 GiB, 4.32% gc time)

@time H = convert_ht(Float32, reffile, trans=true); # 6.059539 seconds (72.59 M allocations: 6.348 GiB, 20.53% gc time)
@time X = convert_gt(Float32, tgtfile, trans=true); # 9.631285 seconds (144.39 M allocations: 12.666 GiB, 18.29% gc time)
@time hapset = compute_optimal_halotype_set(X, H, width = width); # 11.368949 seconds (3.23 M allocations: 1.821 GiB, 3.89% gc time)
@time ph = phase(X, H, hapset=hapset, width=width); # 58.875168 seconds (223.08 k allocations: 36.773 MiB, 0.01% gc time)
@time write_test(ph, H) # 19.748824 seconds (108.38 M allocations: 12.022 GiB, 8.37% gc time)
@time write_test2(ph, H) # 19.439188 seconds (108.38 M allocations: 12.022 GiB, 8.27% gc time)

function write_test(ph, H)
    reader = VCF.Reader(openvcf(tgtfile, "r"))
    writer = VCF.Writer(openvcf(outfile, "w"), header(reader))

    # loop over each record
    for (i, record) in enumerate(reader)
        gtkey = VCF.findgenokey(record, "GT")
        if !isnothing(gtkey)
            # loop over samples
            for (j, geno) in enumerate(record.genotype)
                # if missing = '.' = 0x2e
                if record.data[geno[gtkey][1]] == 0x2e
                    #find where snp is located in phase
                    hap1_position = searchsortedlast(ph[j].strand1.start, i)
                    hap2_position = searchsortedlast(ph[j].strand2.start, i)

                    #find the correct haplotypes
                    hap1 = ph[j].strand1.haplotypelabel[hap1_position]
                    hap2 = ph[j].strand2.haplotypelabel[hap2_position]

                    # save actual allele to data. "0" (REF) => 0x30, "1" (ALT) => 0x31
                    a1, a2 = convert(Bool, H[i, hap1]), convert(Bool, H[i, hap2])
                    record.data[geno[gtkey][1]] = ifelse(a1, 0x31, 0x30)
                    record.data[geno[gtkey][2]] = 0x7c # phased data has separator '|'
                    record.data[geno[gtkey][3]] = ifelse(a2, 0x31, 0x30)
                end
            end
        end
        write(writer, record)
    end

    # close
    flush(writer); close(reader); close(writer)
end

function write_test2(ph, H)
    reader = VCF.Reader(openvcf(tgtfile, "r"))
    writer = VCF.Writer(openvcf(outfile, "w"), header(reader))

    # loop over each record
    for (i, record) in enumerate(reader)
        gtkey = VCF.findgenokey(record, "GT")
        if !isnothing(gtkey)
            # loop over samples
            for (j, geno) in enumerate(record.genotype)
                # if missing = '.' = 0x2e
                if record.data[geno[gtkey][1]] == 0x2e
                    record.data[geno[gtkey][1]] = 0x31
                    record.data[geno[gtkey][2]] = 0x7c # phased data has separator '|'
                    record.data[geno[gtkey][3]] = 0x31
                end
            end
        end
        write(writer, record)
    end

    # close
    flush(writer); close(reader); close(writer)
end



# phase genotype data with beagle 4.1
java -Xss5m -Xmx10g -jar beagle.27Jan18.7e1.jar gt=../AFRped_tgtunphased.vcf ref=../AFRped_ref.vcf niterations=0 out=AFRped_tgtphased





using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using Profile

cd("/Users/biona001/.julia/dev/MendelImpute/simulation")
tgtfile = "target_masked.vcf"
unphase(tgtfile)
X = convert_gt(Float64, tgtfile)
X_unphase = convert_gt(Float64, "unphase." * tgtfile)
all(skipmissing(X .== X_unphase))

compress_vcf_to_gz(tgtfile)
Xgz = convert_gt(Float64, tgtfile * ".gz")
all(skipmissing(X .== Xgz))

compress_vcf_to_gz("unphase.target_masked.vcf")
compress_vcf_to_gz("unphase.target_masked.vcf")



X = convert_gt(Float64, "unphase.target_masked.vcf.gz")
X2 = convert_gt(Float64, "unphase2.target_masked.vcf.gz")




# test if "GT" is allocating a ton
function test()
    s = 0
    for i in 1:100000
        s += get("GT")
    end
    return s
end

get(s) = (s == "GT" ? true : false)

@time test()







# speedup search for best happair

using Revise
using BenchmarkTools
using MendelImpute
using Random
using LinearAlgebra
using LoopVectorization

Random.seed!(123)
n = 500   # number of individuals
p = 400   # number of SNPs in current window
d = 4000  # number of (unique) reference haplotypes
H = convert(Matrix{Float32}, rand(0:1, p, d))
X = convert(Matrix{Float32}, rand(0:2, p, n))
M = Transpose(H) * H
for j in 1:d, i in 1:(j - 1) # off-diagonal
    M[i, j] = 2M[i, j] + M[i, i] + M[j, j]
end
for j in 1:d # diagonal
    M[j, j] *= 4
end
N = Transpose(X) * H
for I in eachindex(N)
    N[I] *= 2
end
happairs = [Tuple{Int, Int}[] for i in 1:n]
hapscore = zeros(eltype(N), n)

function haplopair_test!(
    happairs::Vector{Vector{Tuple{Int, Int}}},
    hapmin::Vector,
    M::AbstractMatrix{T},
    N::AbstractMatrix{T},
    interval::T = convert(T, 3)
    ) where T <: Real

    n, d = size(N)
    fill!(hapmin, typemax(eltype(hapmin)))
    empty!.(happairs)

    @inbounds for k in 1:d, j in 1:k
        # loop over individuals
        @simd for i in 1:n
            score = M[j, k] - N[i, j] - N[i, k]
            # keep all previous best pairs
            if score < hapmin[i]
                push!(happairs[i], (j, k))
                hapmin[i] = score
            end
        end
    end

    return nothing
end

Threads.nthreads()
@time haplopair_test!(happairs, hapscore, M, N)
happairs


# single thread original code
@benchmark haplopair_test!(happairs, hapscore, M, N) # 4.894 s, 4.06 KiB, 1 alloc

# 4 threads, parallelizing by sample
@benchmark haplopair_test!(happairs, hapscore, M, N) # 3.891 s, 6.92 KiB, 31 alloc





using Revise
using GeneticVariation
using Random
using VCFTools
using BenchmarkTools
using MendelImpute

cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf"
marker_chrom, marker_pos, marker_ID, marker_REF, marker_ALT = extract_marker_info(vcffile)







# learn to read/write vcf files faster
using Revise
using GeneticVariation
using Random
using VCFTools
using BenchmarkTools
using MendelImpute

cd("/Users/biona001/.julia/dev/MendelImpute/data/1000_genome_phase3_v5/filtered")

function get_ref_pos(reffile)
    ref_reader = VCF.Reader(openvcf(reffile, "r"))

    ref_marker_pos = zeros(Int, nrecords(reffile))
    for (i, record) in enumerate(ref_reader)
        ref_marker_pos[i] = VCF.pos(record)
    end
    close(ref_reader)

    return ref_marker_pos
end

reffile = "ref.chr20.aligned.vcf.gz" # 379432 entries
@time get_ref_pos(reffile) # 205.682744 seconds (2.75 G allocations: 261.134 GiB, 7.98% gc time)


cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf" # 1356 entries
@time get_ref_pos(vcffile) # 0.686711 seconds (2.39 M allocations: 187.674 MiB, 4.36% gc time)
@time X, pos = convert_gt(Float64, vcffile, trans=true, save_pos=true)



using Revise
using GeneticVariation
using Random
using VCFTools
using BenchmarkTools
using MendelImpute

cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf" # 1356 entries



X, pos = convert_gt(Float64, vcffile, trans=true, save_pos=true)
@code_warntype X, pos = convert_gt(Float64, vcffile, trans=true, save_pos=true)

# type instatiliby caused by varying number of returns
@benchmark convert_gt(Float64, vcffile, trans=true, save_pos=true) # 62.234 ms, 104.77 MiB, 1073112 alloc

# no type instability
@benchmark convert_gt(Float64, vcffile, trans=true) # 68.874 ms, 104.77 MiB, 1073111 alloc




# match indices
X_pos = [41, 51, 76]
H_pos = collect(10:100)
XtoH_idx = indexin(X_pos, H_pos)


i = 1
X_pos[i] == H_pos[XtoH_idx[i]]



X_full = Vector{Union{Missing, Int}}(missing, length(H_pos))
copyto!(@view(X_full[X_pos], [1, 2, 3]))







using LoopVectorization
using BenchmarkTools

function sum_skipmissing(x)
    s = zero(T)
    @inbounds @simd for i in eachindex(x)
        if x[i] !== missing
            s += x[i]
        end
    end
    return s
end

function sum_skipmissing_avx(x)
    s = zero(eltype(x))
    @avx for i in eachindex(x)
        if !ismissing(x[i])
            s += x[i]
        end
    end
    return s
end

function sum_skipmissing_avx2(x)
    s = zero(eltype(x))
    @avx for i in eachindex(x)
        s += ifelse(ismissing(x[i]), zero(eltype(x)), x[i])
    end
    return s
end

x = convert(Vector{Union{Float64, Missing}}, rand(1000));
sum_skipmissing_avx2(x)


@btime euclidean_skipmissing(x, y) #2.727 μs
@btime euclidean_skipmissing_avx(x, y)





# learn to write to file fast
using Random
using CodecZlib
using BenchmarkTools

cd("/Users/biona001/.julia/dev/MendelImpute/test")
function write_test(A)
    io = GzipCompressorStream(open("test.txt.gz", "w"))
    for i in 1:size(A, 1)
        write(io, "hi", "\t", "bye", "\t", "benjamin", "\t", "this is false", "\t", "\t.\tPASS\t.\tGT")
        for j in 1:size(A, 2)
            if A[i, j] < 0.8
                write(io, "\t0/0")
            elseif A[i, j] < 0.9
                write(io, "\t1/0")
            elseif A[i, j] < 1.0
                write(io, "\t1/1")
            else
                error("hi")
            end
        end
        write(io, "\n")
    end
    close(io)
end

function write_test2(A)
    io = GzipCompressorStream(open("test.txt.gz", "w"))
    pb = PipeBuffer()
    for i in 1:size(A, 1)
        print(pb, "hi", "\t", "bye", "\t", "benjamin", "\t", "this is false", "\t", "\t.\tPASS\t.\tGT")
        for j in 1:size(A, 2)
            if A[i, j] < 0.8
                print(pb, "\t0/0")
            elseif A[i, j] < 0.9
                print(pb, "\t1/0")
            elseif A[i, j] < 1.0
                print(pb, "\t1/1")
            else
                error("hi")
            end
        end
        (bytesavailable(pb) > (16*1024)) && write(io, take!(pb))
    end
    write(io, take!(pb))
    close(pb); close(io)
end

Random.seed!(2020)
A = rand(100000, 1000);
# write_test(A)
# write_test2(A)
@time write_test(A)
@time write_test2(A)









function test_lock(b)
    a = 0
    Threads.@threads for i in eachindex(b)
        a += b[i]
    end
    a
end

function test_lock2(b)
    a = 0
    mutex = Threads.SpinLock()
    Threads.@threads for i in eachindex(b)
        lock(mutex)
        a += b[i]
        unlock(mutex)
    end
    a
end

b = 1:1000
test_lock(b)  # 93885
test_lock2(b) # 500500 = correct answer


function text_indexin(a::Vector{Int}, b::Vector{Int}, c::Vector{Int})
    b_in_a = zeros(Int, length(a))
    c_in_a = zeros(Int, length(a))
    for i in 1:2
        b_in_a .= indexin(a, b)
        c_in_a .= indexin(a, c)
    end
    return b_in_a, c_in_a
end





# binary search for vector of invervals
using Revise
using Random
using MendelImpute
using BenchmarkTools

v = [1:4, 5:10, 11:20, 21:35]
searchsortedlast(first.(v), 3)







# generate happairs in windows
using Revise
using Random
using MendelImpute

windows = 10
haplotype_set = [Tuple{Int, Int}[] for i in 1:windows]

Random.seed!(2020)
for w in 1:windows
    haplotype_set[w] = [(rand(1:5), rand(1:5)) for i in 1:rand(1:10)]
end
haplotype_set[2] = haplotype_set[3] = haplotype_set[4] = haplotype_set[5] = haplotype_set[10] = [(-1, -1)]
haplotype_set

originally_empty = falses(windows)
@inbounds for w in 1:windows
    if haplotype_set[w][1] == (-1, -1)
        originally_empty[w] = true
    end
end
originally_empty
nearest_window_with_sufficient_typed_snps(4, haplotype_set, originally_empty)

connect_happairs(haplotype_set)
haplotype_set





(v1.2) pkg> update Distributions
  Updating registry at `~/.julia/registries/General`
  Updating git-repo `https://github.com/JuliaRegistries/General.git`
 Resolving package versions...
 Installed Distributions ─ v0.21.9
  Updating `~/.julia/environments/v1.2/Project.toml`
  [31c24e10] ↑ Distributions v0.21.5 ⇒ v0.21.9
  Updating `~/.julia/environments/v1.2/Manifest.toml`
  [9e28174c] ↑ BinDeps v0.8.10 ⇒ v1.0.1
  [b552c78f] ↑ DiffRules v0.0.10 ⇒ v1.0.1
  [31c24e10] ↑ Distributions v0.21.5 ⇒ v0.21.9
  [276daf66] ↑ SpecialFunctions v0.7.2 ⇒ v0.8.0

(v1.2) pkg> update SpecialFunctions
  Updating registry at `~/.julia/registries/General`
  Updating git-repo `https://github.com/JuliaRegistries/General.git`
ERROR: The following package names could not be resolved:
 * SpecialFunctions (276daf66-3868-5448-9aa4-cd146d93841b in manifest but not in project)
Please specify by known `name=uuid`.

(v1.2) pkg> update VCFTools SnpArrays
  Updating registry at `~/.julia/registries/General`
  Updating git-repo `https://github.com/JuliaRegistries/General.git`
  Updating git-repo `https://github.com/OpenMendel/VCFTools.jl`
  Updating git-repo `https://github.com/OpenMendel/SnpArrays.jl`
 Resolving package versions...
  Updating `~/.julia/environments/v1.2/Project.toml`
  [4e780e97] ~ SnpArrays v0.3.2 #master (https://github.com/OpenMendel/SnpArrays.jl)
  [a620830f] ~ VCFTools v0.1.2 #master (https://github.com/OpenMendel/VCFTools.jl)
  Updating `~/.julia/environments/v1.2/Manifest.toml`
  [4e780e97] ~ SnpArrays v0.3.2 #master (https://github.com/OpenMendel/SnpArrays.jl)
  [a620830f] ~ VCFTools v0.1.2 #master (https://github.com/OpenMendel/VCFTools.jl)

julia>
(base) [biona001@n2238 ~]$ ls
haplotype_comparisons  IHT  ProPCA  R  thyrosim  yeast_eQTL  yeast_eQTL_newton
(base) [biona001@n2238 ~]$ cd haplotype_comparisons/HRC/
(base) [biona001@n2238 HRC]$ clear


using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using StatsBase
using CodecZlib
using ProgressMeter
using JLD2, FileIO, JLSO
using BenchmarkTools
using GroupSlices

cd("/u/home/b/biona001/haplotype_comparisons/HRC")
function filter_and_mask(chr::Int, maf::Float64)
    # filter chromosome data for unique snps
    println("filtering for unique snps")
    data = "HRC.r1-1.EGA.GRCh37.chr$chr.haplotypes.vcf.gz"
    full_record_index = .!find_duplicate_marker(data)
    @time VCFTools.filter(data, full_record_index, 1:nsamples(data),
        des = "chr$chr.uniqueSNPs.vcf.gz")

    # summarize data
    println("summarizing data")
    total_snps, samples, _, _, _, maf_by_record, _ = gtstats("chr$chr.uniqueSNPs.vcf.gz")

    # generate target panel with all snps
    println("generating complete target panel")
    n = 1000
    sample_idx = falses(samples)
    sample_idx[1:n] .= true
    shuffle!(sample_idx)
    @time VCFTools.filter("chr$chr.uniqueSNPs.vcf.gz", 1:total_snps,
        sample_idx, des = "target.chr$chr.full.vcf.gz", allow_multiallelic=false)

    # also generate reference panel without target samples
    println("generating reference panel without target samples")
    @time VCFTools.filter("chr$chr.uniqueSNPs.vcf.gz", 1:total_snps,
        .!sample_idx, des = "ref.chr$chr.excludeTarget.vcf.gz", allow_multiallelic=false)

    # generate target file with 1000 samples and typed snps with certain maf
    println("generating target file with typed snps only")
    my_maf = findall(x -> x > maf, maf_by_record)
    p = length(my_maf)
    record_idx = falses(total_snps)
    record_idx[my_maf] .= true
    @time VCFTools.filter("chr$chr.uniqueSNPs.vcf.gz", record_idx, sample_idx,
        des = "target.chr$chr.typedOnly.maf$maf.vcf.gz", allow_multiallelic=false)

    # unphase and mask 1% entries in target file
    println("unphasing and masking entries in target file with typed snps only")
    masks = falses(p, n)
    missingprop = 0.01
    for j in 1:n, i in 1:p
        rand() < missingprop && (masks[i, j] = true)
    end
    @time mask_gt("target.chr$chr.typedOnly.maf$maf.vcf.gz", masks,
        des="target.chr$chr.typedOnly.maf$maf.masked.vcf.gz", unphase=true)

    # finally compress reference file to jlso format
    reffile = "ref.chr$chr.excludeTarget.vcf.gz"
    outfile = "ref.chr$chr.excludeTarget.jlso"
    @time compress_haplotypes(reffile, outfile, width=2048);
end
Random.seed!(2020)
chr = ARGS[1]
maf = ARGS[2]
println("Running chrom $chr. Typed SNPs have frequency maf > $maf")
fdsa
@time filter_and_mask(chr, maf)




using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using StatsBase
using CodecZlib
using ProgressMeter
using JLD2, FileIO, JLSO
using BenchmarkTools
using GroupSlices

chr = parse(Int, ARGS[1])
maf = parse(Float64, ARGS[2])
width = 2048

# searching single bkpts, deleting suboptimal pairs in dp, skip windows with <50 typed snps
Random.seed!(2020)
tgtfile = "target.chr$chr.typedOnly.maf$maf.masked.vcf.gz"
reffile = "ref.chr$chr.excludeTarget.jlso"
outfile = "mendel.imputed.dp$width.chr$chr.maf$maf.vcf.gz"
@time ph, hs = phase(tgtfile, reffile, outfile=outfile, impute=true, width=width)

X_complete = convert_gt(UInt8, "target.chr$chr.full.vcf.gz")
n, p = size(X_complete)
X_mendel = convert_gt(UInt8, outfile)
println("error overall = $(sum(X_mendel .!= X_complete) / n / p) \n")




# try different window width to see number of haplotypes per window

using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using StatsBase
using CodecZlib
using ProgressMeter
using JLD2, FileIO, JLSO
using BenchmarkTools
using GroupSlices

vcffile = "ref.chr21.excludeTarget.vcf.gz"
widths  = [128, 256, 512, 1024]
dims = 2
flankwidth = 0
trans = (dims == 2 ? true : error("currently VCFTools only import chr/pos...info when H transposed"))
H, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt = convert_ht(Bool, vcffile, trans=trans, save_snp_info=true, msg="importing vcf data...")

for width in widths
    outfile = "ref.chr21.w$width.excludeTarget.jlso"

    snps = (dims == 2 ? size(H, 1) : size(H, 2))
    windows = floor(Int, snps / width)
    compressed_Hunique = MendelImpute.CompressedHaplotypes(windows, width, snps, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt)

    avg_uniqhaps = 0
    for w in 1:windows
        if w == 1
            cur_range = 1:(width + flankwidth)
        elseif w == windows
            cur_range = ((windows - 1) * width - flankwidth + 1):snps
        else
            cur_range = ((w - 1) * width - flankwidth + 1):(w * width + flankwidth)
        end
        compressed_Hunique.CWrange[w] = cur_range

        H_cur_window = (dims == 2 ? view(H, cur_range, :) : view(H, :, cur_range))
        hapmap = groupslices(H_cur_window, dims=dims)
        unique_idx = unique(hapmap)
        complete_to_unique = indexin(hapmap, unique_idx)
        uniqueH = (dims == 2 ? H_cur_window[:, unique_idx] : H_cur_window[unique_idx, :])
        compressed_Hunique[w] = MendelImpute.CompressedWindow(unique_idx, hapmap, complete_to_unique, uniqueH)
        avg_uniqhaps += length(unique_idx)
    end
    avg_uniqhaps /= windows
    println("Width = $width averages $avg_uniqhaps haplotype per window. Total win = $windows")

    # save using JLSO or JLD2
    endswith(outfile, ".jld2") && JLD2.@save outfile compressed_Hunique
    endswith(outfile, ".jlso") && JLSO.save(outfile, :compressed_Hunique => compressed_Hunique, format=:julia_serialize, compression=:gzip)

end




function test(storage, a, b)
    UpperTriangular(storage) .= UpperTriangular(a .* b')
end

storage = zeros(1000, 1000);
a = rand(1000);
b = rand(1000);

@time test(storage, a, b);




using BenchmarkTools, Random, Lasso

"""
Find minimum of `M_jk - N_ij - N_ik` for each i. M is upper triangular.

Saves minimum value in `best_err` and `j, k` index in `row` and `col`.
"""
function best_index!(
    best_err::Vector{T},
    row::Vector{Int32},
    col::Vector{Int32},
    M::AbstractMatrix{T},
    N::AbstractMatrix{T},
    ) where T <: Real

    n, d = size(N)

    @inbounds for k in 1:d, j in 1:k
        @simd for i in 1:n
            score = M[j, k] - N[i, j] - N[i, k]

            # save minimum value and index
            if score < best_err[i]
                best_err[i], row[i], col[i] = score, j, k
            end
        end
    end
end

Random.seed!(2020)
n = 100
p = 100
d = 5000

row_index, col_index = ones(Int32, n), ones(Int32, n)
best_err = [typemax(Float32) for _ in 1:n]
M = rand(Float32, d, d)
N = rand(Float32, n, d)

# use Int64: 886.923 ms (0 allocations: 0 bytes)
@btime best_index!($best_err, $row_index, $col_index, $M, $N) setup=(fill!(best_err, Inf32); row_index = ones(Int, n); col_index = ones(Int, n))

# use Int32: 892.310 ms (0 allocations: 0 bytes)
@btime best_index!($best_err, $row_index, $col_index, $M, $N) setup=(fill!(best_err, Inf32); row_index = ones(Int32, n); col_index = ones(Int32, n))


using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using StatsBase
using CodecZlib
using ProgressMeter
using JLD2, FileIO, JLSO
using BenchmarkTools
using GroupSlices

tgtfile = "target.chr18.typedOnly.maf0.1.masked.vcf.gz"
reffile = "ref.chr18.excludeTarget.vcf.gz"
outfile = "ref.chr18.excludeTarget.jlso"
width = 512

compress_haplotypes(reffile, tgtfile, outfile, width)




using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using StatsBase
using CodecZlib
using ProgressMeter
using JLD2, FileIO, JLSO
using BenchmarkTools
using GroupSlices

Random.seed!(2020)
n = 1000
p = 512
d = 1000

X = rand(UInt8, p, n)
H = bitrand(p, d)

@btime haplopair_thin(X, H, keep=100)  #224.270 ms (5478 allocations: 8.17 MiB)
@btime haplopair_thin2(X, H, keep=100) #69.389 ms (5018 allocations: 15.55 MiB)


@btime haplopair_thin(X, H, keep=600)  #1.708 s (6507 allocations: 10.50 MiB)
@btime haplopair_thin2(X, H, keep=600) #1.159 s (5019 allocations: 15.55 MiB)






using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using StatsBase
using CodecZlib
using ProgressMeter
using JLD2, FileIO, JLSO
using BenchmarkTools
using GroupSlices
using Lasso

Random.seed!(2020)
n = 1000
p = 512
d = 1000

X = rand(p, n)
H = rand(0.:1., p, d)
happairs, hapscore = haplopair(X, H)

# try lasso on x1
x1 = X[:, 1];
path = fit(LassoPath, H, x1, Normal(), IdentityLink(), λ=collect(0.035:0.001:0.039))
path.coefs[:, 3] # picks 335 and 443
happairs[1][1], happairs[2][1] # our code picks (453, 543)

# try lasso on x2
x2 = X[:, 2];
path = fit(LassoPath, H, x2, Normal(), IdentityLink(), λ=collect(0.035:0.001:0.039))
path.coefs[:, 3] # picks 139 and 983
happairs[1][2], happairs[2][2] # our code picks (428, 518)


# find λ so that only 2 β is non-zero
function bisection()
    λtop = 1.0
    λlow = 0.0
    while true
        λmid = (λtop + λlow) / 2
        path = fit(LassoPath, H, x2, Normal(), IdentityLink(), λ=[λmid])

        nz = count(!iszero, path.coefs)
        if nz == 2
            println(path.coefs.rowval)
            println("λtop = $λtop, λlow = $λlow, λmid = $λmid, nz = $nz")
            break
        else
            # bisection search
            nz < 2 ? (λtop = λmid) : (λlow = λmid)
        end
        println("λtop = $λtop, λlow = $λlow, λmid = $λmid, nz = $nz")
    end
end
bisection()




using MendelImpute, JLSO, UnicodePlots, StatsBase

@time loaded = JLSO.load("ref.chr20.w64.maf0.01.excludeTarget.jlso") # roughly 30 sec
compressed_Hunique = loaded[:compressed_Hunique];
x = count_haplotypes_per_window(compressed_Hunique)
UnicodePlots.histogram(x, nbins=50)





using Revise
using VCFTools
using MendelImpute
using JLSO, JLD2, FileIO
using LinearAlgebra
cd("/home/biona001/HRC")

# to load:
@load "example_data.jld2" Ms Ns
Ms[1] # M for width = 64
Ms[2] # M for width = 128
Ns[1] # N for width = 64 ... etc


# width = 64 averages 873 haplotype per window. Total win =2789
# width = 128 averages 2497 haplotype per window. Total win = 1394
# width = 256 averages 6985 haplotype per window. Total win =697
# width = 512 averages 17234 haplotype per window. Total win = 348
# width = 1024 averages 32428 haplotype per window. Total win = 174

function create_M_and_N()
    w = 1 # choose a window
    widths = [64, 128, 256, 512, 1024]
    Ms = Vector{Matrix{Float32}}(undef, 5)
    Ns = Vector{Matrix{Float32}}(undef, 5)

    # genotype data
    tgtfile = "target.chr20.typedOnly.maf0.01.masked.vcf.gz"
    X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt = VCFTools.convert_gt(UInt8, tgtfile, trans=true, save_snp_info=true, msg = "Importing genotype file...")

    for (i, width) in enumerate(widths)
        # load data (takes roughly 30 sec)
        @time loaded = JLSO.load("ref.chr20.w$width.maf0.01.excludeTarget.jlso")
        compressed_Hunique = loaded[:compressed_Hunique];
        Hw_aligned = compressed_Hunique.CW_typed[w].uniqueH
        Xw_idx_start = (w - 1) * width + 1
        Xw_idx_end = (w == windows ? length(X_pos) : w * width)
        Xw_aligned = X[Xw_idx_start:Xw_idx_end, :]

        # haplopair and haplopair! in haplotype_pair.jl
        p, n  = size(Xw_aligned)
        d     = size(Hw_aligned, 2)
        Xwork = zeros(Float32, p, n)
        Hwork = convert(Matrix{Float32}, Hw_aligned)
        MendelImpute.initXfloat!(Xw_aligned, Xwork) # initializes missing
        M = zeros(Float32, d, d)
        N = zeros(Float32, n, d)

        # assemble M
        mul!(M, Transpose(Hwork), Hwork)
        for j in 1:d, i in 1:(j - 1) # off-diagonal
            M[i, j] = 2M[i, j] + M[i, i] + M[j, j]
        end
        for j in 1:d # diagonal
            M[j, j] *= 4
        end
        # assemble N
        mul!(N, Transpose(Xwork), Hwork)
        @simd for I in eachindex(N)
            N[I] *= 2
        end

        Ms[i] = M
        Ns[i] = N
    end
    return Ms, Ns
end
Ms, Ns = create_M_and_N()
@save "example_data.jld2" Ms Ns



using Random
using SparseArrays
using DataFrames

Random.seed!(2020)
x = sprandn(1_000_000, 10_000, 0.1)
I, J, V = findnz(x)
df = DataFrame([:I => I, :J => J, :V => V])
CSV.write("/tmp/spmatrix.csv", df)



using Random
using BenchmarkTools

function set_index!(x, idx)
    x .= false
    for i in idx
        x[i] = true
    end
end

x = bitrand(10000)
idx = rand(1:10000, 1000)

@benchmark set_index!($x, $idx) setup=(idx = rand(1:10000, 1000))



struct bitv
    v::BitVector
end

function set_index!(x::bitv, idx)
    x.v .= false
    @inbounds for i in idx
        x.v[i] = true
    end
end
x = bitv(bitrand(10000))
idx = rand(1:10000, 1000)

@benchmark set_index!($x, $idx) setup=(idx = rand(1:10000, 1000))




using Random
using BenchmarkTools

function noelapse(n)
    s = 0.0
    for i in 1:n
        s += rand()
    end
    s
end

function elapse(n)
    s = 0.0
    t = @elapsed begin
        for i in 1:n
            s += rand()
        end
    end
    s, t
end

function elapse2(n)
    s = 0.0
    t = 0.0
    @inbounds for i in 1:n
        t += @elapsed s += rand()
    end
    s, t
end

@benchmark x = noelapse(100)
@benchmark x, t = elapse(100)
@benchmark x, t = elapse2(100)







using Random
using BenchmarkTools

# see post https://discourse.julialang.org/t/editing-an-array-with-multiple-threads/25364/9
function f(spacing)
    n = 1_000_000
    x = rand(n)

    s = zeros(Threads.nthreads()*spacing) # sum of x
    c = zeros(Threads.nthreads()*spacing) # count additions in each thread

    Threads.@threads for i = 1:n
        s[Threads.threadid()*spacing] += x[i] # thread i adds to s[i]
        c[Threads.threadid()*spacing] += 1    # thread i adds to c[i]
    end
    return sum(s), c
end

@btime f(1); #17.915 ms (35 allocations: 7.63 MiB)
@btime f(8); #1.513 ms (35 allocations: 7.63 MiB)



# Phase genotypes with beagle 5.1
java -jar ../beagle.18May20.d20.jar gt=target.typedOnly.maf0.01.masked.vcf.gz ref=ref.excludeTarget.vcf.gz out=beagle.target.phased nthreads=8 impute=false

# run tabix to get tabix index for .vcf.gz files
tabix -p vcf beagle.target.phased.vcf.gz
tabix -p vcf ref.excludeTarget.vcf.gz

# Run impute 5
impute5 --h ref.excludeTarget.vcf.gz --m ../maps/genetic_maps.b37.tar/chr20.b37.gmap.gz --g beagle.target.phased.vcf.gz --r 20 --o impute5.imputed.vcf.gz


let 
    tot = 0
    for (key, val) in compressed_Hunique.CW_typed[1].hapmap
        if length(val) > 196
            tot += 1
        end
    end
    println(tot)
end


using Revise
using LinearAlgebra
using BenchmarkTools
using ThreadPools

BLAS.set_num_threads(1)

function test_multiply()
    threads = Threads.nthreads()
    A = [rand(1000, 1000) for _ in 1:threads]
    B = [rand(1000, 1000) for _ in 1:threads]
    C = [rand(1001, 1001) for _ in 1:threads]
    # Threads.@threads for i in 1:10000
    ThreadPools.@qthreads for i in 1:10000
            id = Threads.threadid()
        mul!(view(C[id], 1:1000, 1:1000), A[id], B[id])
    end
    return C
end
test_multiply();






using Revise
using VCFTools
using GroupSlices

cd(".julia/dev/MendelImpute/simulation/compare2")
H = convert_ht(Bool, "ref.excludeTarget.vcf.gz", trans=true)
result = get_window_intervals(H, 1000)


# experiment with parallel write
using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using StatsBase
using CodecZlib
using ProgressMeter
using BenchmarkTools
using GroupSlices

function Base.write(
    outfile::AbstractString,
    X::AbstractMatrix,
    )
    # write minimal meta information to outfile
    io = openvcf(outfile, "w")
    pb = PipeBuffer()
    print(pb, "##fileformat=VCFv4.2\n")
    print(pb, "##source=MendelImpute\n")
    print(pb, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

    pmeter = Progress(size(X, 1), 5, "Writing to file...")
    @inbounds for i in 1:size(X, 1)
        # print ith record
        write_snp!(pb, @view(X[i, :]))

        (bytesavailable(pb) > (1024*1024)) && write(io, take!(pb))
        next!(pmeter)
    end
    write(io, take!(pb))

    # close & return
    close(io); close(pb)
    return nothing
end
"""
Helper function for saving a record (SNP), not tracking phase information.
"""
function write_snp!(pb::IOBuffer, X::AbstractVector)
    n = length(X)
    @inbounds for j in 1:n
        if X[j] == 0
            print(pb, "\t0/0")
        elseif X[j] == 1
            print(pb, "\t1/0")
        elseif X[j] == 2
            print(pb, "\t1/1")
        else
            error("imputed genotypes can only be 0, 1, 2 but got $(X[j])")
        end
    end
    print(pb, "\n")
    nothing
end

function write_threaded(
    outfile::AbstractString,
    X::AbstractMatrix,
    )
    threads = Threads.nthreads()
    snps = size(X, 1)
    len = div(snps, threads)
    files = ["tmp$i.vcf.gz" for i in 1:threads]

    # write minimal meta information to outfile
    io = [openvcf(files[i], "w") for i in 1:threads]
    pb = [PipeBuffer() for _ in 1:threads]
    print(pb[1], "##fileformat=VCFv4.2\n")
    print(pb[1], "##source=MendelImpute\n")
    print(pb[1], "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    pmeter = Progress(snps, 5, "Writing to file...")

    # each thread writes `len` SNPs
    Threads.@threads for t in 1:threads
        id = Threads.threadid()
        cur_ranges = (id == threads ? ((threads-1)*len+1:snps) : (1:len) .+ (t-1)*len)
        @inbounds for i in cur_ranges
            write_snp!(pb[id], @view(X[i, :]))
            (bytesavailable(pb[id]) > (1024*1024)) && write(io[id], take!(pb[id])) 
            next!(pmeter)
        end
        write(io[id], take!(pb[id]))
    end
    close.(io); close.(pb) # close io and buffer

    # concatenate all files into 1 VCF file
    run(pipeline(`cat $files`, stdout=outfile))

    # delete intermediate files
    for i in 1:threads
        rm("tmp$i.vcf.gz", force=true)
    end

    return nothing
end


Random.seed!(2020)
p = 100_000
n = 1000
X = convert(Matrix{UInt8}, rand(0:2, p, n))
Y = convert(Matrix{UInt8}, rand(0:2, 10, 10))
outfile = "test.vcf.gz"

@time write(outfile, X) # 23.534699 seconds (60.28 k allocations: 384.933 MiB, 0.39% gc time)
@time write_threaded(outfile, X) #


# minimal example
using CodecZlib
for i in 1:3
    file = "file$i.gz"
    x = rand(2, 2)
    io = GzipCompressorStream(open(file, "w"))
    write(io, "hello, testing from i = $i")
    # close(io)
end

files = ["file1.gz", "file2.gz", "file3.gz"]
run(pipeline(`cat $files`, stdout="allfiles.gz"))




using ProfileView
@profview write(outfile, Y) 
@profview write(outfile, X) 


# experiment with parallel read
using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using StatsBase
using CodecZlib
using ProgressMeter
using BenchmarkTools
using GroupSlices


struct test
    x::Int
    y::Int
end
hi = test(1, 2)
@code_warntype hi.x

struct test2
    x::Int
    y::Int
    z::Vector{Int}
end
hi = test2(1, 2, collect(1:3))
@code_warntype hi.x

struct test3
    x::Int
    y::Float64
end
hi = test3(1, 2.0)
@code_warntype hi.x

struct test5{T, U}
    x::T
    y::U
end
hi = test5(1, 2.0)
@code_warntype hi.x

struct CW
    uniqueindex::Vector{Int32}
    hapmap::Dict{Int32, Vector{Int32}}
    to_unique::Vector{Int32}
    uniqueH::BitMatrix
end
fd = CW(rand(Int32, 0:4))



# import per SNP quality score (postimputation)
using GeneticVariation
using VCFTools

cd("/Users/biona001/.julia/dev/MendelImpute/simulation/compare2")
vcf = "mendel.imputed.vcf.gz"
reader = VCF.Reader(openvcf(vcf, "r"))




# test BGEN
using BGEN
b = Bgen(BGEN.datadir("example.8bits.bgen"); 
    sample_path=BGEN.datadir("example.sample"), 
    idx_path=BGEN.datadir("example.8bits.bgen.bgi"))
b = Bgen(BGEN.datadir("haplotypes.bgen"); 
    idx_path=BGEN.datadir("haplotypes.bgen.bgi"))

variants = parse_variants(b; from_bgen_start=true);
minor_allele_dosage!(b, variants[1]; T=Float32)
probabilities!(b, variants[1]; T=Float32)
v = variants[1]

function convert_gt(b::Bgen, T=Float32)
    n = n_samples(b)
    p = n_variants(b)

    # return arrays
    G = Matrix{T}(undef, p, n)
    GsampleID = Vector{String}(undef, n)
    Gchr = Vector{String}(undef, p)
    Gpos = Vector{String}(undef, p)
    GsnpID = Vector{String}(undef, p)
    Gref = Vector{String}(undef, p)
    Galt = Vector{String}(undef, p)

    # loop over each variant
    i = 1
    for v in iterator(b; from_bgen_start=true)
        dose = minor_allele_dosage!(b, v; T=T)
        copyto!(@view(G[i, :]), dose)
        chr[i], pos[i], snpID[i], ref[i], alt[i] = chrom(v), pos(v), rsid(v),
            major_allele(v), minor_allele(v)
        i += 1
        clear!(v)
    end
    return G, sampleID, chr, pos, snpID, ref, alt
end
G, sampleID, chr, pos, snpID, ref, alt = convert_gt(b)



b2 = Bgen(BGEN.datadir("complex.bgen"))
vs = parse_variants(b2)
p = probabilities!(b2, vs[1]) #ploidy 1

n_samples(b2)
n_variants(b2)


b = Bgen(BGEN.datadir("example.8bits.bgen"); 
    sample_path=BGEN.datadir("example.sample"), 
    idx_path=BGEN.datadir("example.8bits.bgen.bgi"))



using Revise, MendelImpute
using BGEN, VCFTools, BGEN
cd("/Users/biona001/Desktop/qctool_v2.0.8-osx")

nsamples("ref.excludeTarget.vcf"), nrecords("ref.excludeTarget.vcf")
b = Bgen("ref.excludeTarget.bgen")
variants = parse_variants(b; from_bgen_start=true);
v = variants[1]
p = probabilities!(b, v)
n = 2n_samples(b)
p = n_variants(b)
H = BitMatrix(undef, p, n)
Hchr = Vector{String}(undef, p)
Hpos = Vector{Int}(undef, p)
HsnpID = Vector{String}(undef, p)
Href = Vector{String}(undef, p)
Halt = Vector{String}(undef, p)


Gtest = MendelImpute.import_target("ref.excludeTarget.bgen")

Htest = MendelImpute.import_reference("ref.excludeTarget.bgen")
Htrue = convert_ht(Bool, "ref.excludeTarget.vcf", trans=true)