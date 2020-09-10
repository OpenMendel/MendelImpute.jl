@testset "search breakpoints" begin
    Random.seed!(2019)
    p, d = 1000, 20
    H = rand(0.0:1.0, p, d)

    bkpt = p >>> 1 # 500
    X = Vector{Union{Float64, Missing}}(undef, p)
    mask = [rand() < 0.1 for i in 1:p] 
    # 1 | 2
    # 1 | 3
    X .= H[:, 1] + [H[1:bkpt, 2]; H[bkpt+1:p, 3]]
    X[mask] .= missing
    bkpt_optim, err_optim = MendelImpute.search_breakpoint(X, H[:, 1], H[:, 2], H[:, 3])
    # @code_warntype MendelImpute.search_breakpoint(X, H[:, 1], H[:, 2], H[:, 3])
    @test bkpt_optim == 499
    @test err_optim == 0

    # 4 | 5
    # 2 | 5
    bkpt = 200
    X .= H[:, 5] + [H[1:bkpt, 4]; H[bkpt+1:p, 2]]
    X[mask] .= missing
    bkpt_optim, err_optim = MendelImpute.search_breakpoint(X, H[:, 5], H[:, 4], H[:, 2])
    @test bkpt_optim == 198
    @test err_optim == 0

    # 2 | 3
    # 5 | 4
    bkpt1 = 200
    bkpt2 = 500
    X .= [H[1:bkpt1, 2]; H[bkpt1+1:p, 5]] + [H[1:bkpt2, 3]; H[bkpt2+1:p, 4]]
    X[mask] .= missing
    bkpt_optim, err_optim = MendelImpute.search_breakpoint(X, H[:, 2], H[:, 5], H[:, 3], H[:, 4])
    @test bkpt_optim == (199, 499)
    @test err_optim == 0

    # 5 | 3
    # 15 | 14
    bkpt1 = 300
    bkpt2 = 700
    X .= [H[1:bkpt1, 5]; H[bkpt1+1:p, 15]] + [H[1:bkpt2, 3]; H[bkpt2+1:p, 14]]
    X[mask] .= missing
    bkpt_optim, err_optim = MendelImpute.search_breakpoint(X, H[:, 5], H[:, 15], H[:, 3], H[:, 14])
    @test bkpt_optim == (297, 699)
    @test err_optim == 0
end