module HaplotypingTest

using MendelImpute
using Base.Test
using BenchmarkTools

n, p, d = 1000, 100, 20
H = rand(0:1, d, p)
X = rand(0:2, n, p)
#@code_warntype haplopair(X, H)
n, d = size(X, 1), size(H, 1)
M = A_mul_Bt(H, H)
M .*= 2
N = A_mul_Bt(X, H)
N .*= 2
happair  = zeros(Int, n, 2)
hapscore = zeros(eltype(N), n)
@time haplopair!(happair, hapscore, M, N)

end
