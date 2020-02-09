###### GOAL: #######
### By exchanging bits of 2 bitvectors at the same position, 
### minimize the number of breakpoints in both strands, and
### ruturn that number.
####################
# Example 1
# x = [0 1 1 1]
# y = [1 0 0 0]
# We can exchange position 1:
# x = [1 1 1 1]
# y = [0 0 0 0]
# then breakpoint = 0
# 
# Example 2
# x = [0 1 0 0 1]
# y = [1 0 1 0 0]
# we can exchange position 2:
# x = [0 0 0 0 1]
# y = [1 1 1 0 0]
# then breakpoint = 2 (1 at end of x and another at middle of y)

using Random
using Test
using BenchmarkTools

"""
Starting at position `n`, search through to the end of `seq1` and `seq2` to see if 
position `seq1[n - 1]` and `seq2[n - 1]` should be swapped. Swapping index are 
stored in `flip_idx`. 

Returns the total number of breakpoints in `seq1` and `seq2` after swapping. 
"""
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
        # only flip bits in previous position
        return binary_flip!(n + 1, seq1, seq2, intermediates)
    elseif seq1[n] == seq2[n]
        # don't flip and only calculate error
        return binary_flip!(n + 1, seq1, seq2, flip_idx, intermediates) + !(seq1[n] == seq1[n - 1]) + !(seq2[n] == seq2[n - 1])
    else
        # calculate error of flip/noflip by recursion
        yesflip = binary_flip!(n + 1, seq1, seq2, flip_idx, intermediates) + !(seq1[n] == seq2[n - 1]) + !(seq2[n] == seq1[n - 1])
        noflip  = binary_flip!(n + 1, seq1, seq2, flip_idx, intermediates) + !(seq1[n] == seq1[n - 1]) + !(seq2[n] == seq2[n - 1])
        # println("n = $n, yesflip = $yesflip, noflip = $noflip, min = $(min(yesflip, noflip))")
        # store intermediate results for later retrival
        intermediates[n] = min(yesflip, noflip)
        if yesflip < noflip
            # record flipping location, flip 2 sequence at previous location, 
            flip_idx[n - 1] = true
            seq1[n - 1], seq2[n - 1] = seq2[n - 1], seq1[n - 1]
            # println("flipped $(n - 1)!")
        end
        return min(yesflip, noflip)
    end
end
binary_flip!(seq1, seq2, flip) = binary_flip!(2, seq1, seq2, flip) #start at position 2


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

# nonmemoized (bitrand(100) would fail)
x = bitrand(30)
y = bitrand(30)
flip = falses(30)
@btime binary_flip!($x, $y, $flip) # 298.914 μs (1 allocation: 336 bytes)

# memoized 
x = bitrand(30)
y = bitrand(30)
@btime binary_flip!($x, $y, $flip) # 737.567 ns (1 allocation: 336 bytes)

x = bitrand(10000)
y = bitrand(10000)
@time binary_flip!(x, y)


Random.seed!(2020)
x = bitrand(10)
y = bitrand(10)
[x y]
flip = falses(10);
@time binary_flip!(x, y, flip)
[x y]

Random.seed!(2020)
x = bitrand(20)
y = bitrand(20)
[x y]
flip = falses(20);
@time binary_flip!(x, y, flip)
[x y]

#more examples
x = bitrand(10)
y = bitrand(10)
[x y]
flip = falses(10);
@time binary_flip!(x, y, flip)
[x y]

x = bitrand(10000)
y = bitrand(10000)
[x y]
flip = falses(10000);
@time binary_flip!(x, y, flip)
[x y]

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

function has_intersect!(c::BitVector, a::BitVector, b::BitVector)
    c .= a .& b
    return any(c)
end

function has_intersect(a::BitVector, b::BitVector)
    @inbounds for i in eachindex(a)
        if a[i] && b[i]
            return true
        end
    end
    return false
end

a = falses(100_000)
b = falses(100_000)
c = falses(100_000)

Random.seed!(2020)
a[1:1000] .= b[1:1000] .= true
shuffle!(a)
shuffle!(b)

@btime has_intersect($a, $b)
@btime has_intersect!($c, $a, $b)

function set_flip!(
    n::Int, 
    vector_set1::Vector{BitVector}, 
    vector_set2::Vector{BitVector}, 
    flip_idx::BitVector,
    intermediates::Vector{Int} = [-1 for i in 1:(length(vector_set1) + 1)],
    storage::BitVector = falses(length(vector_set1[1]))
    )
    # quick lookup
    if intermediates[n] != -1
        return intermediates[n]
    end

    if n > length(vector_set1)
        return 0
    elseif n == 1
        # only flip bits in previous position
        return set_flip!(n + 1, vector_set1, vector_set2, flip_idx, intermediates)
    elseif vector_set1[n] == vector_set2[n]
        # don't flip and only calculate error
        return set_flip!(n + 1, vector_set1, vector_set2, flip_idx, intermediates) + 
                !has_intersect!(storage, vector_set1[n], vector_set1[n - 1]) +
                !has_intersect!(storage, vector_set2[n], vector_set2[n - 1])
    else
        # calculate error of flip/noflip by recursion
        yesflip = set_flip!(n + 1, vector_set1, vector_set2, flip_idx, intermediates) + 
                    !has_intersect!(storage, vector_set1[n], vector_set2[n - 1]) + 
                    !has_intersect!(storage, vector_set2[n], vector_set1[n - 1])
        noflip  = set_flip!(n + 1, vector_set1, vector_set2, flip_idx, intermediates) + 
                    !has_intersect!(storage, vector_set1[n], vector_set1[n - 1]) + 
                    !has_intersect!(storage, vector_set2[n], vector_set2[n - 1])

        # store intermediate results for later retrival
        intermediates[n] = min(yesflip, noflip)
        if yesflip < noflip
            # record flipping location, flip 2 sequence at previous location, 
            flip_idx[n - 1] = true
            vector_set1[n - 1], vector_set2[n - 1] = vector_set2[n - 1], vector_set1[n - 1]
        end
        return min(yesflip, noflip)
    end
end
set_flip!(vector_set1, vector_set2, flip) = set_flip!(2, vector_set1, vector_set2, flip) #start at position 2


function has_intersect!(c::BitVector, a::BitVector, b::BitVector)
    c .= a .& b
    return any(c)
end

using Random
using Test
using BenchmarkTools

Random.seed!(2020)

x = [bitrand(10), bitrand(10), bitrand(10), bitrand(10)]
y = [bitrand(10), bitrand(10), bitrand(10), bitrand(10)]
flip = falses(4);
@test set_flip!(x, y, flip) == 2
@test flip == [false; true; false; false]


x = [bitrand(5), bitrand(5), bitrand(5), bitrand(5), bitrand(5), bitrand(5), bitrand(5)]
y = [bitrand(5), bitrand(5), bitrand(5), bitrand(5), bitrand(5), bitrand(5), bitrand(5)]
flip = falses(7);
@test set_flip!(x, y, flip) == 1
@test flip == [false; false; false; false; true; false; false]







using Random
using Test
using BenchmarkTools

Random.seed!(1)

a = bitrand(5)
b = bitrand(5)

struct bitvec_pair

end

for pair in (a, b)
    println(pair)
end




