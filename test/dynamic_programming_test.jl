@testset "binary_flip!" begin
    x = [1 0]
    y = [0 1]
    flip = falses(2)
    @test binary_flip!(x, y, flip) == 0
    @test flip == [true; false]
    @test x == [0 0]
    @test y == [1 1]

    x = [1 0 1]
    y = [0 1 0]
    flip = falses(3)
    @test binary_flip!(x, y, flip) == 0
    @test flip == [false; true; false]
    @test x == [1 1 1]
    @test y == [0 0 0]

    x = [0 1 1 1]
    y = [1 0 0 0]
    flip = falses(4)
    @test binary_flip!(x, y, flip) == 0
    @test flip == [true; false; false; false]
    @test x == [1 1 1 1]
    @test y == [0 0 0 0]

    x = [1 0 1 0]
    y = [0 1 0 0]
    flip = falses(4)
    @test binary_flip!(x, y, flip) == 1
    @test flip == [false; true; false; false]
    @test x == [1 1 1 0]
    @test y == [0 0 0 0]

    x = [0 1 0 0 1]
    y = [1 0 1 0 0]
    flip = falses(5)
    @test binary_flip!(x, y, flip) == 2
    @test flip == [false; true; false; false; false]
    @test x == [0 0 0 0 1]
    @test y == [1 1 1 0 0]

    x = [1 0 1 0]
    y = [0 1 0 1]
    flip = falses(4)
    @test binary_flip!(x, y, flip) == 0
    @test flip == [true; false; true; false]
    @test x == [0 0 0 0]
    @test y == [1 1 1 1]

    x = [0 1 1 1 0 0]
    y = [1 1 1 0 0 0]
    flip = falses(6)
    @test binary_flip!(x, y, flip) == 3
    @test flip == [false; false; false; false; false; false]
    @test x == [0 1 1 1 0 0]
    @test y == [1 1 1 0 0 0]

    x = [0 1 0 1 1 0]
    y = [1 0 1 0 0 0]
    flip = falses(6)
    @test binary_flip!(x, y, flip) == 1
    @test flip == [true; false; true; false; false; false]
    @test x == [1 1 1 1 1 0]
    @test y == [0 0 0 0 0 0]

    x = [0 1 1 1 1 0]
    y = [1 0 1 0 0 0]
    flip = falses(6)
    @test binary_flip!(x, y, flip) == 3
    @test flip == [true; false; false; false; false; false]
    @test x == [1 1 1 1 1 0]
    @test y == [0 0 1 0 0 0]

    x = [0 1 1 1 1 0]
    y = [1 0 1 0 0 1]
    flip = falses(6)
    @test binary_flip!(x, y, flip) == 2
    @test flip == [true; false; false; true; true; false]
    @test x == [1 1 1 0 0 0]
    @test y == [0 0 1 1 1 1]

    x = [0 1 1 1 1 0]
    y = [1 0 1 1 0 1]
    flip = falses(6)
    @test binary_flip!(x, y, flip) == 2
    @test flip == [true; false; false; false; true; false]
    @test x == [1 1 1 1 0 0]
    @test y == [0 0 1 1 1 1]
end

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

    @test memory[3][(1, 1)] == 2
    @test memory[3][(1, 5)] == 1
    @test memory[3][(3, 1)] == 2
    @test memory[3][(3, 5)] == 1

    @test memory[2][(1, 5)] == 1
    @test memory[2][(1, 7)] == 2
    @test memory[2][(1, 8)] == 2
    @test memory[2][(2, 5)] == 2
    @test memory[2][(2, 7)] == 3
    @test memory[2][(2, 8)] == 3
    @test memory[2][(6, 5)] == 2
    @test memory[2][(6, 7)] == 3
    @test memory[2][(6, 8)] == 3

    @test memory[1][(1, 4)] == 2
    @test memory[1][(1, 5)] == 1
    @test memory[1][(1, 6)] == 2
    @test memory[1][(2, 4)] == 3
    @test memory[1][(2, 5)] == 2
    @test memory[1][(2, 6)] == 3
    @test memory[1][(3, 4)] == 3
    @test memory[1][(3, 5)] == 2
    @test memory[1][(3, 6)] == 3

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

    @test memory[4][(5, 3)] == 1
    @test memory[4][(5, 9)] == 1
    @test memory[4][(3, 9)] == 0
    @test memory[4][(8, 10)] == 1

    @test memory[3][(1, 5)] == 2
    @test memory[3][(6, 9)] == 1

    @test memory[2][(2, 9)] == 2
    @test memory[2][(3, 4)] == 3
    @test memory[2][(1, 2)] == 3
    @test memory[2][(8, 8)] == 3
    @test memory[2][(6, 1)] == 2

    @test memory[1][(3, 10)] == 4
end