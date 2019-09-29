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

testset "redundant_haplotypes" begin
	#import data
	cd("/Users/biona001/.julia/dev/MendelImpute/test")
	rawdata = readdlm("AFRped_geno.txt", ',', Float32);
	people = 664;
	X = copy(Transpose(rawdata[1:people, 1:(end - 1)]));
	function create_hap(x)
	    n, p = size(x)
	    h = one(eltype(x))
	    for j in 1:p, i in 1:n
	        if x[i, j] != 0
	            x[i, j] -= h
	        end
	    end
	    return copy(Transpose(x))
	end
	H = create_hap(rawdata[(people + 1):end, 1:(end - 1)]);

	#mask random entries
	Random.seed!(123)
	missingprop = 0.1
	p, n = size(X)
	X2 = Matrix{Union{Missing, eltype(X)}}(X)
	Xm = ifelse.(rand(eltype(X), p, n) .< missingprop, missing, X2)
	Xm_original = copy(Xm)
	width = 64
	windows = floor(Int, p / width)

	hapset = redundant_haplotypes(Xm, H, width=width)

	#check if 10th window is correct
	w = 10
    cur_range = ((w - 1) * width + 1):(w * width)
	H_cur = H[cur_range, :]	
	result = collect(hapset.strand1.p[w, 1]) 
	all_col_should_agree = H_cur[:, result]
	@test all(all_col_should_agree[:, 1] .== all_col_should_agree[:, 2])
	@test all(all_col_should_agree[:, 4] .== all_col_should_agree[:, 8])
	@test all(all_col_should_agree[:, 2] .== all_col_should_agree[:, 7])

	# check if last window is correct
	last_win_length = mod(p, width)
	Hlast = H[(end - last_win_length):end, :]
	result = collect(hapset.strand1.p[end, 1])
	all_col_should_agree = Hlast[:, result]
	@test all(all_col_should_agree[:, 1] .== all_col_should_agree[:, 8])
	@test all(all_col_should_agree[:, 3] .== all_col_should_agree[:, 4])
	@test all(all_col_should_agree[:, 2] .== all_col_should_agree[:, 7])
 end