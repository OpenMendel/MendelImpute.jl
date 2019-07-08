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
@benchmark unique(H, dims=1)               #
@benchmark filter_redundant_haplotypes(H)  #






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




