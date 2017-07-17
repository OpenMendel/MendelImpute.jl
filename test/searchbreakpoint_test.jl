module HaplotypingTest

using MendelImpute, NullableArrays
using Base.Test, BenchmarkTools

srand(123)
p, d = 1000, 100
H = rand(0.:1., p, d)
# 1 | 2
# 1 | 3
bkpt = round(Int, p / 2)
X = NullableArray(H[:, 1] + [H[1:bkpt, 2]; H[bkpt+1:p, 3]], bitrand(p))

info("search optimal break point in 1 strand")
#@code_warntype search_breakpoint(X, H, 1, (2, 3))
bkpt_optim, err_optim = search_breakpoint(X, H, 1, (2, 3))
@show bkpt_optim, err_optim
@time search_breakpoint(X, H, 1, (2, 3))
@time search_breakpoint(X, H, 1, (2, 3))

#@code_warntype search_breakpoint_new(X, H, 1, (2, 3))
bkpt_optim, err_optim = search_breakpoint_new(X, H, 1, (2, 3))
@show bkpt_optim, err_optim
@time search_breakpoint_new(X, H, 1, (2, 3))

# info("search optimal break point in 2 strands")
# #@code_warntype search_breakpoint(X, H, (1, 2), (3, 4))
# bkpts_optim, err_optim = search_breakpoint(X, H, (1, 4), (2, 3))
# for bkpt1 in 0:10
#     for bkpt2 in 0:10
#         errors = countnz([H[1:bkpt1, 1]; H[bkpt1+1:10, 4]] +
#             [H[1:bkpt2, 2]; H[bkpt2+1:10, 3]] .≠ X)
#         println("optimal breakpoints ($bkpt1, $bkpt2), $errors errors")
#     end
# end
# @show bkpts_optim, err_optim
#
# info("continue haplotypes")
# #@code_warntype continue_haplotype(X, H, happair_prev, happair_next, bkpt)
# # both strands match
# happair_next_optim, breakpts = continue_haplotype(X, H, (1, 2), (1, 2))
# @test happair_next_optim == (1, 2)
# @test breakpts == (-1, -1)
# happair_next_optim, breakpts = continue_haplotype(X, H, (1, 2), (2, 1))
# @test happair_next_optim == (1, 2)
# @test breakpts == (-1, -1)
# # one strand matches
# happair_next_optim, breakpts = continue_haplotype(X, H, (1, 2), (1, 3))
# @test happair_next_optim == (1, 3)
# happair_next_optim, breakpts = continue_haplotype(X, H, (1, 2), (3, 1))
# @test happair_next_optim == (1, 3)
# happair_next_optim, breakpts = continue_haplotype(X, H, (1, 2), (3, 2))
# @test happair_next_optim == (3, 2)
# happair_next_optim, breakpts = continue_haplotype(X, H, (1, 2), (2, 3))
# @test happair_next_optim == (3, 2)
# # no strand matches
# happair_next_optim, breakpts = continue_haplotype(X, H, (1, 2), (3, 4))
# @show happair_next_optim, breakpts

end
