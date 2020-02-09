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