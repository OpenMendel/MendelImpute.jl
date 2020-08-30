@testset "search breakpoints" begin
	Random.seed!(2019)
	p, d = 1000, 20
	H = rand(0.:1., p, d)
	bkpt = p >>> 1 # 500
	X = Vector{Union{Float64, Missing}}(undef, p)
	mask = [rand() < 0.1 for i in 1:p] 
	# 1 | 2
	# 1 | 3
	X .= H[:, 1] + [H[1:bkpt, 2]; H[bkpt+1:p, 3]]
	X[mask] .= missing
	bkpt_optim, err_optim = MendelImpute.search_breakpoint(X, H[:, 1], H[:, 2], H[:, 3])

	# @code_warntype search_breakpoint(X, H[:, 1], H[:, 2], H[:, 3])
	@test bkpt_optim == 499
	@test err_optim == 0
 
	# 2 | 3
	# 2 | 4
	bkpt = 800
	X .= H[:, 2] + [H[1:bkpt, 3]; H[bkpt+1:p, 4]]
	X[mask] .= missing
	bkpt_optim, err_optim = MendelImpute.search_breakpoint(X, H[:, 2], H[:, 3], H[:, 4])

	@test bkpt_optim == 800
	@test err_optim == 0

	# 4 | 5
	# 2 | 5
	bkpt = 200
	X .= H[:, 5] + [H[1:bkpt, 4]; H[bkpt+1:p, 2]]
	X[mask] .= missing
	bkpt_optim, err_optim = MendelImpute.search_breakpoint(X, H[:, 5], H[:, 4], H[:, 2])

	@test bkpt_optim == 198
	@test err_optim == 0
end