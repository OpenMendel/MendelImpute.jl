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

function binary_flip(n, seq1, seq2, 
    intermediates::Vector{Int} = [-1 for i in 1:(length(seq1) + 1)]
    )
    # quick lookup
    if intermediates[n] != -1
        return intermediates[n]
    end

    if n > length(seq1)
        return 0
    elseif n == 1
        return binary_flip(n + 1, seq1, seq2, intermediates)
    elseif seq1[n] == seq2[n]
        return !(seq1[n] == seq1[n - 1]) + !(seq2[n] == seq2[n - 1]) + binary_flip(n + 1, seq1, seq2, intermediates)
    else
        yesflip = !(seq1[n] == seq2[n - 1]) + !(seq2[n] == seq1[n - 1]) + binary_flip(n + 1, seq1, seq2, intermediates)
        noflip  = !(seq1[n] == seq1[n - 1]) + !(seq2[n] == seq2[n - 1]) + binary_flip(n + 1, seq1, seq2, intermediates)
        # println("n = $n, yesflip = $yesflip, noflip = $noflip, min = $(min(yesflip, noflip))")
        intermediates[n] = min(yesflip, noflip)
        return min(yesflip, noflip)
    end
end
binary_flip(seq1, seq2) = binary_flip(2, seq1, seq2) #start at position 2

x = [0 1 1 1]
y = [1 0 0 0]
@test binary_flip(x, y) == 0

x = [0 1 1 0]
y = [1 0 0 0]
@test binary_flip(x, y) == 1

x = [0 1 1 1 0 0]
y = [1 1 1 0 0 0]
@test binary_flip(x, y) == 3

x = [0 1 0 1 1 0]
y = [1 0 1 0 0 0]
@test binary_flip(x, y) == 1

x = [0 1 1 1 1 0]
y = [1 0 1 0 0 0]
@test binary_flip(x, y) == 3

# nonmemoized (bitrand(100) would fail)
x = bitrand(30)
y = bitrand(30)
@btime binary_flip($x, $y) # 54.562 μs (1 allocation: 336 bytes)

# memoized 
x = bitrand(30)
y = bitrand(30)
@btime binary_flip($x, $y) # 484.790 ns (1 allocation: 336 bytes)

