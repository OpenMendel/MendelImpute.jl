testset "unique_haplotypes" begin
    H = zeros(Int, 1, 9)
	H[1] = 1
	H[2] = 2
	H[3] = 2
	H[4] = 3
	H[5] = 2
	H[6] = 4
	H[7] = 1
	H[8] = 4
	H[9] = 1
	hapset = unique_haplotypes(H, 128, 'T')

	@test all(hapset.uniqH[1]  .== [1, 2, 4, 6])
	@test all(hapset.hapmap[1] .== [1, 2, 2, 4, 2, 6, 1, 6, 1])
end