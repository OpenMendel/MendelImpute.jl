@testset "connect_happairs" begin
    T = Tuple{Int, Int}

    # first case
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

    sol_path, memory, best_err = connect_happairs(haplotype_set)

    @test best_err == 1.0
    @test sol_path == [(1, 5), (1, 5), (1, 5), (5, 4)]

    @test memory[3][(1, 1)] == (2.0, (4, 4))
    @test memory[3][(1, 5)] == (1.0, (5, 4))
    @test memory[3][(3, 1)] == (2.0, (4, 4))
    @test memory[3][(3, 5)] == (1.0, (5, 4))

    @test memory[2][(1, 5)] == (1.0, (1, 5))
    @test memory[2][(1, 7)] == (2.0, (1, 5))
    @test memory[2][(1, 8)] == (2.0, (1, 5))
    @test memory[2][(2, 5)] == (2.0, (1, 5))
    @test memory[2][(2, 7)] == (3.0, (1, 5))
    @test memory[2][(2, 8)] == (3.0, (1, 5))
    @test memory[2][(6, 5)] == (2.0, (1, 5))
    @test memory[2][(6, 7)] == (3.0, (1, 5))
    @test memory[2][(6, 8)] == (3.0, (1, 5))

    @test memory[1][(1, 4)] == (2.0, (1, 5))
    @test memory[1][(1, 5)] == (1.0, (1, 5))
    @test memory[1][(1, 6)] == (2.0, (1, 5))
    @test memory[1][(2, 4)] == (3.0, (1, 5))
    @test memory[1][(2, 5)] == (2.0, (1, 5))
    @test memory[1][(2, 6)] == (3.0, (1, 5))
    @test memory[1][(3, 4)] == (3.0, (1, 5))
    @test memory[1][(3, 5)] == (2.0, (1, 5))
    @test memory[1][(3, 6)] == (3.0, (1, 5))

    # second case
    windows = 5
    haplotype_set = [T[] for i in 1:windows]

    Random.seed!(2020)
    for w in 1:windows
        haplotype_set[w] = [(rand(1:10), rand(1:10)) for i in 1:rand(1:10)]
    end
    sol_path, memory, best_err = connect_happairs(haplotype_set)

    @test best_err == 4.0
    @test sol_path == [(3, 10), (2, 9), (6, 9), (3, 9), (9, 3)]

    @test memory[4][(5, 3)] == (1.0, (9, 3))
    @test memory[4][(5, 9)] == (1.0, (9, 3))
    @test memory[4][(3, 9)] == (0.0, (9, 3))
    @test memory[4][(8, 10)] == (1.0, (1, 10))

    @test memory[3][(1, 5)] == (2.0, (5, 3))
    @test memory[3][(6, 9)] == (1.0, (3, 9))

    @test memory[2][(2, 9)] == (2.0, (6, 9))
    @test memory[2][(3, 4)] == (3.0, (6, 9))
    @test memory[2][(1, 2)] == (3.0, (1, 5))
    @test memory[2][(8, 8)] == (3.0, (6, 9))
    @test memory[2][(6, 1)] == (2.0, (6, 9))

    @test memory[1][(3, 10)] == (4.0, (2, 9))
end