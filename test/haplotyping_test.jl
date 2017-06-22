module HaplotypingTest

using MendelImpute
using Base.Test

info("test isuniquerows")

srand(123)
A = rand(0:1, 8, 4)
isuniq = zeros(Bool, size(A, 1))
p = collect(1:size(A, 1))
finduniquerows!(isuniq, p, A)
@show A
@show isuniq
@show p
@show view(A, isuniq, :)

Av = view(A, 1:2:7, 1:2:3)
finduniquerows!(zeros(Bool, size(Av, 1)), collect(1:size(Av, 1)), Av)


end
