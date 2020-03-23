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

    sol_path, next_pair, subtree_err, best_err = connect_happairs(haplotype_set)

    @test best_err == 1.0
    @test sol_path == [(1, 5), (1, 5), (1, 5), (5, 4)]
    @test all(subtree_err[1] .== [2.0; 1.0; 2.0; 3.0; 2.0; 3.0; 3.0; 2.0; 3.0])
    @test all(subtree_err[2] .== [1.0; 2.0; 2.0; 2.0; 3.0; 3.0; 2.0; 3.0; 3.0])
    @test all(subtree_err[3] .== [2.0; 1.0; 2.0; 1.0])
    @test all(subtree_err[4] .== [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0])
    @test all(next_pair[1] .== [1; 1; 1; 1; 1; 1; 1; 1; 1])
    @test all(next_pair[2] .== [2; 2; 2; 2; 2; 2; 2; 2; 2])
    @test all(next_pair[3] .== [1; 4; 1; 4])
    @test all(next_pair[4] .== [0; 0; 0; 0; 0; 0; 0; 0; 0])

    # second case
    windows = 5
    haplotype_set = [T[] for i in 1:windows]

    Random.seed!(2020)
    for w in 1:windows
        haplotype_set[w] = [(rand(1:10), rand(1:10)) for i in 1:rand(1:10)]
    end
    sol_path, next_pair, subtree_err, best_err = connect_happairs(haplotype_set)

    @test best_err == 4.0
    @test sol_path == [(3, 10), (2, 9), (6, 9), (3, 9), (9, 3)]
    @test all(subtree_err[1] .== [4.0])
    @test all(subtree_err[2] .== [2.0; 3.0; 3.0; 2.0; 3.0; 2.0])
    @test all(subtree_err[3] .== [2.0; 1.0])
    @test all(subtree_err[4] .== [1.0; 1.0; 0.0; 1.0])
    @test all(subtree_err[5] .== [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0;])
    @test all(next_pair[1] .== [1])
    @test all(next_pair[2] .== [2; 2; 1; 2; 2; 2])
    @test all(next_pair[3] .== [1; 3])
    @test all(next_pair[4] .== [1; 1; 1; 3])
    @test all(next_pair[5] .== [0; 0; 0; 0; 0; 0; 0])
end