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
        write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        for i in 1:Int(d / 2)
            write(io, "ref1$i\t")
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
        write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        for i in 1:d
            write(io, "target$i\t")
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
    simulate_correlated_snparray(n, p, s; block_length, hap, prob)

Simulates a SnpArray with correlation. SNPs are divided into blocks where each
adjacent SNP is the same with probability prob. There are no correlation between blocks.

# Arguments:
- `p`: number of SNPs
- `d`: number of people
- `s`: name of SnpArray that will be created (memory mapped) in the current directory. To not memory map, use `undef`.

# Output:
* `H`: `p x d` haplotype matrix. Each column is a haplotype.

# Optional arguments:
- `block_length`: length of each LD block
- `hap`: haplotype pool size for each block
- `prob`: with probability `prob` an adjacent SNP would be the same. 
"""
function simulate_haplotypes(p::Int64, d::Int64, s::Union{String, UndefInitializer}; 
            block_length::Int64=20, hap::Int=20, prob::Float64=0.75)
    
    @assert mod(p, block_length) == 0 "block_length ($block_length) is not divible by p ($p)"
    @assert 0 < prob < 1 "transition probably should be between 0 and 1, got $prob"

    H = BitArray(undef, p, d)
    haplotypes = zeros(hap, block_length)
    snps = zeros(block_length)
    blocks = Int(p / block_length)

    @inbounds for b in 1:blocks

        #create pool of haplotypes for each block
        _sample_haptotypes!(haplotypes, prob)

        for i in 1:n
            #sample 2 haplotypes with replacement from the pool of haplotypes
            row1 = rand(1:hap)
            row2 = rand(1:hap)
            for j in 1:block_length
                snps[j] = haplotypes[row1, j] + haplotypes[row2, j]
            end

            #copy haplotypes into x
            _copy_blocks!(x, i, snps, b, block_length)
        end
    end

    return x
end

function _sample_haptotypes!(haplotypes::Matrix, prob::Float64)
    n, p = size(haplotypes)
    fill!(haplotypes, 0)

    @inbounds for i in 1:n
        cur_row_sum = 0
        while cur_row_sum == 0
            curr = rand(0:1)
            haplotypes[i, 1] = curr
            cur_row_sum += curr
            for j in 2:p
                stay = rand(Bernoulli(prob)) #stay = 1 means retain the current value
                curr = (stay == 1 ? curr : 1 - curr)
                haplotypes[i, j] = curr
                cur_row_sum += curr
            end
        end
    end
end

# function _copy_blocks!(x::SnpArray, row, snps, cur_block, block_length)
#     #copy sampled snps into SnpArray
#     @inbounds for k in 1:length(snps)
#         c = snps[k]
#         col = (cur_block - 1) * block_length + k
#         if c == 0
#             x[row, col] = 0x00
#         elseif c == 1
#             x[row, col] = 0x02
#         elseif c == 2
#             x[row, col] = 0x03
#         else
#             throw(error("SNP values should be 0, 1, or 2 but was $c"))
#         end
#     end
# end