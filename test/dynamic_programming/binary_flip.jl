###### GOAL: #######
### Determine minimum number of breakpoints in both strands.
### You can exchange bits in 2 bitvectors at the same position.
####################
# e.g. 
# x = [0 1 1 1]
# y = [1 0 0 0]
# breakpoint(x, y) = 2 (since there is flip from position 1 to 2 in both)
# But we can exchange position 1:
# x = [1 1 1 1]
# y = [0 0 0 0]
# then breakpoint(x, y) = 0

using Random
using Test
using BenchmarkTools

function binary_flip!(n, seq1, seq2, flip_idx,
    intermediates::Vector{Int} = [-1 for i in 1:(length(seq1) + 1)]
    )
    # quick lookup
    if intermediates[n] != -1
        return intermediates[n]
    end

    if n > length(seq1)
        return 0
    elseif n == 1
        return binary_flip!(n + 1, seq1, seq2, intermediates)
    elseif seq1[n] == seq2[n]
        return binary_flip!(n + 1, seq1, seq2, flip_idx, intermediates) + !(seq1[n] == seq1[n - 1]) + !(seq2[n] == seq2[n - 1])
    else
        yesflip = binary_flip!(n + 1, seq1, seq2, flip_idx, intermediates) + !(seq1[n] == seq2[n - 1]) + !(seq2[n] == seq1[n - 1])
        noflip  = binary_flip!(n + 1, seq1, seq2, flip_idx, intermediates) + !(seq1[n] == seq1[n - 1]) + !(seq2[n] == seq2[n - 1])
        # println("n = $n, yesflip = $yesflip, noflip = $noflip, min = $(min(yesflip, noflip))")
        intermediates[n] = min(yesflip, noflip)
        if yesflip < noflip
            flip_idx[n - 1] = true
            seq1[n - 1], seq2[n - 1] = seq2[n - 1], seq1[n - 1]
            yesflip, noflip = noflip, yesflip
            println("flipped $(n - 1)!")
        end
        return min(yesflip, noflip)
    end
end
binary_flip!(seq1, seq2, flip) = binary_flip!(2, seq1, seq2, flip) #start at position 2

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
@test x = [1 1 1 1]
@test y = [0 0 0 0]

x = [1 0 1 0]
y = [0 1 0 0]
flip = falses(4)
@test binary_flip!(x, y, flip) == 1
@test flip == [false; true; false; false]
@test x == [1 1 1 0]
@test y == [0 0 0 0]

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

# nonmemoized (bitrand(100) would fail)
x = bitrand(30)
y = bitrand(30)
@btime binary_flip($x, $y) # 54.562 μs (1 allocation: 336 bytes)

# memoized 
x = bitrand(30)
y = bitrand(30)
@btime binary_flip($x, $y) # 484.790 ns (1 allocation: 336 bytes)

x = bitrand(10000)
y = bitrand(10000)
@time binary_flip(x, y)






function a()
    println("reached a!")
    return 1
end
function b()
    println("reached b!")
    return 2
end
c = a() + b() 
# right hand side of equations are evaluated from left to right
# julia> c = a() + b()
# reached a!
# reached b!
# 3
