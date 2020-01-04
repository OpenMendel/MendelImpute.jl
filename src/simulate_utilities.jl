"""
    make_vcf_file(H, filename, phased)

Creates a .vcf file with `filename` based on reference panels `H`. Consecutive columns are treated as genotypes

# Inputs:
+ `H`: Matrix of 0s or 1s. Each column is a haplotype. e.g. Columns 1 and 2 form the genotype for person 1. 
+ `filename`: A string for the resulting .vcf file. 
+ `phased`: True uses '|' as separator. False uses '/' as separator. 
"""
function make_refvcf_file(
    H::BitArray{2};
    filename="simulated_ref.vcf", 
    phased = true
    )

    p, d = size(H)
    separator = (phased ? '|' : '/')
    iseven(d) || error("make_vcf_file: number of haplotypes must be even but was $d")
    endswith(filename, ".vcf") || (filename = filename * ".vcf")

    open(filename, "w") do io
        # minimal meta information
        write(io, "##fileformat=VCFv4.3\n")
        write(io, "##source=MendelImpute\n")
        write(io, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

        # header line
        write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for i in 1:Int(d / 2)
            write(io, "\tref$i")
        end
        write(io, "\n")

        # write phase info
        for snp in 1:p
            write(io, "1\t$snp\tref_snp_$snp\tA\tC\t.\tPASS\t.\tGT")
            for i in 1:2:d
                write(io, string("\t", Int(H[snp, i]), separator, Int(H[snp, i + 1])))
            end
            write(io, "\n")
        end
    end
end

function make_tgtvcf_file(
    X::Matrix{Union{Int, Missing}};
    filename="simulated_tgt.vcf", 
    )

    p, d = size(X)
    endswith(filename, ".vcf") || (filename = filename * ".vcf")

    open(filename, "w") do io
        # minimal meta information
        write(io, "##fileformat=VCFv4.3\n")
        write(io, "##source=MendelImpute\n")
        write(io, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

        # header line
        write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for i in 1:d
            write(io, "\ttarget$i")
        end
        write(io, "\n")

        # write phase info
        for snp in 1:p
            write(io, "1\t$snp\ttgt_snp_$snp\tA\tC\t.\tPASS\t.\tGT")
            for i in 1:d
                if ismissing(X[snp, i])
                    genotype = "./."
                elseif X[snp, i] == 2
                    genotype = "1/1"
                elseif X[snp, i] == 1
                    genotype = "1/0"
                elseif X[snp, i] == 0
                    genotype = "0/0"
                else
                    error("genotypes can only be 0, 1, 2, or missing, but was $(X[snp, i])")
                end
                write(io, string("\t", genotype))
            end
            write(io, "\n")
        end
    end
end

"""
    simulate_markov_haplotypes(p, d, prob)

Simulates a haplotype matrix as a markov chain. The `i`th allele (0 or 1) will transition 
to the opposite allele (0 or 1) with probability `prob` at the `i + 1`th allele. 

# Inputs
- `p`: Length of each haplotype.
- `d`: Total number of haplotypes. 
- `prob`: transition probability.

# Output
- `H`: `p x d` haplotype matrix. Each column is a haplotype
"""
function simulate_markov_haplotypes(
    p::Int64, 
    d::Int64;
    prob = 0.25,
    )
    @assert 0 < prob < 1 "transition probably `prob` should be between 0 and 1, got $prob"

    H = falses(p, d)
    @inbounds for j in 1:d
        H[1, j] = rand(Bool)
        for i in 2:p
            H[i, j] = (rand() < prob ? !H[i - 1, j] : H[i - 1, j])
        end
    end
    return H
end

"""
    simulate_uniform_haplotypes(p, d, prob)

Simulates a haplotype matrix `H` where `H[i, j] = 1` with probability `prob`. 

# Inputs
- `p`: Length of each haplotype.
- `d`: Total number of haplotypes. 
- `prob`: probability that an entry in H is 1.

# Output
- `H`: `p x d` haplotype matrix. Each column is a haplotype
"""
function simulate_uniform_haplotypes(
    p::Int64, 
    d::Int64;
    prob = 0.25,
    )
    @assert 0 < prob < 1 "prob should be between 0 and 1, got $prob"

    H = falses(p, d)
    @inbounds for j in 1:d, i in 1:p
        if rand() < prob
            H[i, j] = true
        end
    end
    return H
end

"""
    simulate_genotypes(H; block_length)

Simulates a genotype matrix `X` from a pool of haplotypes `H`. SNPs are 
divided into blocks of length `block_length` where a pair of haplotypes 
are drawn uniformly from `H` to form the genotype of that given block. 

# Arguments:
- `H`: `p x d` haplotype matrix. Each column is a haplotype. 
- `people`: number of samples desired for the genotype matrix
- `block_length`: length of each LD block

# Output:
* `X`: `p x people` genotype matrix. Each column is a person's genotype. 
"""
function simulate_genotypes(
    H::BitArray{2}; 
    people::Int = size(H, 2),
    block_length::Int64=113
    )
    
    p, d = size(H)
    X = zeros(Int, p, people)
    blocks = Int(ceil(p / block_length))

    # for each block, sample 2 ` with replacement from the pool of haplotypes
    for b in 1:(blocks - 1), i in 1:people
        hap1 = rand(1:d)
        hap2 = rand(1:d)
        block_start = (b - 1) * block_length
        for j in 1:block_length
            X[block_start + j, i] = H[block_start + j, hap1] + H[block_start + j, hap2]
        end
    end

    # treat last block separately
    for i in 1:people
        hap1 = rand(1:d)
        hap2 = rand(1:d)
        block_start = (blocks - 1) * block_length
        for j in 1:(p - block_start)
            X[block_start + j, i] = H[block_start + j, hap1] + H[block_start + j, hap2]
        end
    end

    return X
end
