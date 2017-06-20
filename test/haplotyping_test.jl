module HaplotypingTest

using MendelImpute
using Base.Test

info("test isuniquerows")

srand(123)
A = rand(0:1, 8, 4)
isuniq = falses(size(A, 1))
p = ones(Int, size(A, 1))
isuniquerows!(isuniq, p, A)
@show A
@show isuniq
@show p
@show view(A, isuniq, :)

end
