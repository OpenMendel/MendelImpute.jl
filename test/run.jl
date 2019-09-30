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
H = bitrand(1280, 1000)
# Hwork = unique_haplotypes(H, 1:128)
# Hunique = unique_haplotypes(H, 128)

for i in 1:500
	H[:, 2i] .= H[:, 2i - 1]
end

# Hwork = unique_haplotypes(H, 1:128, 'T')
Hunique = unique_haplotypes(H, 128, 'T')

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
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using MendelImpute
using Random
using Profile
using ElasticArrays

cd("/Users/biona001/.julia/dev/MendelImpute/test")

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
width = 400
windows = floor(Int, p / width)

copyto!(Xm, Xm_original)

#Hua's code error = 0.013899062884015356 (width = 64)
#Hua's code error = 0.0033518189729195617 (width = 400) #12.170407 seconds (87.59 k allocations: 23.879 MiB, 0.15% gc time)
hapset = phase(Xm, H, width)
impute2!(Xm, H, hapset) 

#current code error = 0.01967560196500439 (width = 64)
#current code error = 0.009237186022036839 (width = 3*64)
#current code error = 0.006945615284063669 (width = 400)
#current code error = 0.004905503108811067 (width = 3*400) #3.769867 seconds (663.78 k allocations: 61.506 MiB, 0.80% gc time)
copyto!(Xm, Xm_original)
hapset = phase2(Xm, H, width=width)

copyto!(Xm, Xm_original)
hapset = phase2(Xm, H, width=3*width)

#non-redundant haplotypes without search breakpoints has error = 0.08077261261736622
hapset = phase3(Xm, H, width=width)


# look at the haplotype intersections
hapset.strand1.p[1:10, 1]
hapset.strand2.p[1:10, 1]


# calculate error rate
missing_idx    = ismissing.(Xm_original)
total_missing  = sum(missing_idx)
actual_missing_values  = convert(Vector{Int64}, X[missing_idx])  #true values of missing entries
imputed_missing_values = convert(Vector{Int64}, Xm[missing_idx]) #imputed values of missing entries
error_rate = sum(actual_missing_values .!= imputed_missing_values) / total_missing



function naive_impute!(X)
	n, p = size(X)
    fillval = convert(eltype(X), 1.0)
	for j in 1:p, i in 1:n
		ismissing(X[i, j]) && (X[i, j] = fillval)
	end
end
naive_impute!(Xm) #fill with 0 gives error rate = 0.09844272948788609
naive_impute!(Xm) #fill with 1 gives error rate = 0.9343025343313078
naive_impute!(Xm) #fill with 2 gives error rate = 0.9672547361808062
















