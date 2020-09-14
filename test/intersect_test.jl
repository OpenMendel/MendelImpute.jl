@testset "intersect_lange!" begin
    x = [1, 3, 4, 5, 7, 9]
    y = [2, 3, 5, 6]
    @test MendelImpute.intersect_size_sorted(x, y) == 2
    MendelImpute.intersect_lange!(x, y)
    @test all(x .== [3, 5])
    @test all(y .== [2, 3, 5, 6])

    x = [1, 3, 4, 7]
    y = [2, 3, 5, 6, 7, 10]
    @test MendelImpute.intersect_size_sorted(x, y) == 2
    MendelImpute.intersect_lange!(x, y)
    @test all(x .== [3, 7])
    @test all(y .== [2, 3, 5, 6, 7, 10])

    # allow repeats, although we don't have any in MendelImpute
    x = [3, 4, 7, 7, 7, 10] 
    y = [2, 3, 5, 7, 7, 10]
    @test MendelImpute.intersect_size_sorted(x, y) == 4
    MendelImpute.intersect_lange!(x, y)
    @test all(x .== [3, 7, 7, 10])
    @test all(y .== [2, 3, 5, 7, 7, 10])
end