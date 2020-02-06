using Random

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

function binary_flip(n, seq1, seq2)
    length(seq1) == length(seq2) || error("seq1 and seq2 have different length")

    if n == 0
        return 0
    elseif seq1[n] == seq2[n]
        return !(seq1[n] == seq1[n + 1]) + !(seq2[n] == seq2[n + 1]) + binary_flip(n - 1, seq1, seq2)
    else
        yesflip = !(seq1[n] == seq2[n + 1]) + !(seq2[n] == seq1[n + 1]) + binary_flip(n - 1, seq1, seq2)
        noflip  = !(seq1[n] == seq1[n + 1]) + !(seq2[n] == seq2[n + 1]) + binary_flip(n - 1, seq1, seq2)
        # println("n = $n, yesflip = $yesflip, noflip = $noflip, min = $(min(yesflip, noflip))")
        return min(yesflip, noflip)
    end
end

x = [0 1 1 1]
y = [1 0 0 0]
binary_flip(3, x, y) # should be 0

x = [0 1 1 0]
y = [1 0 0 0]
binary_flip(3, x, y) # should be 1

x = [0 1 1 1 0 0]
y = [1 1 1 0 0 0]
binary_flip(5, x, y) # should be 3

x = [0 1 0 1 1 0]
y = [1 0 1 0 0 0]
binary_flip(5, x, y) # should be 1

x = [0 1 1 1 1 0]
y = [1 0 1 0 0 0]
binary_flip(5, x, y) # should be 3



