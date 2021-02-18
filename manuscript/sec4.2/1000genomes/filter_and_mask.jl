using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using JLSO
using LinearAlgebra

function filter_and_mask(data, maf, chr)
    # filter chromosome data for unique snps
    #println("filtering for unique snps"); flush(stdout)
    #full_record_index = .!find_duplicate_marker(data)
    #@time VCFTools.filter(data, full_record_index, 1:nsamples(data), 
    #    des = "chr$chr.uniqueSNPs.vcf.gz")

    # summarize data
    println("summarizing data"); flush(stdout)
    total_snps, samples, _, _, _, maf_by_record, _ = gtstats("chr$chr.uniqueSNPs.vcf.gz")

    # keep snps with at least 5 copies of minor alleles
    snps_tokeep = findall(x -> x ≥ 5 / 2samples, maf_by_record)

    # generate target panel with all snps
    println("generating complete target panel")
    n = 100
    sample_idx = falses(samples)
    sample_idx[1:n] .= true
    shuffle!(sample_idx)
    @time VCFTools.filter("chr$chr.uniqueSNPs.vcf.gz", snps_tokeep, 
        sample_idx, des = "target.chr$chr.full.vcf.gz", allow_multiallelic=false)

    #also generate reference panel without target samples
    println("generating reference panel without target samples")
    @time VCFTools.filter("chr$chr.uniqueSNPs.vcf.gz", snps_tokeep, 
        .!sample_idx, des = "ref.chr$chr.excludeTarget.vcf.gz", allow_multiallelic=false)

    # generate target file with 100 samples and typed snps with certain maf
    println("generating target file with typed snps only")
    my_maf = findall(x -> x > maf, maf_by_record)  
    p = length(my_maf)
    record_idx = falses(total_snps)
    record_idx[my_maf] .= true
    @time VCFTools.filter("chr$chr.uniqueSNPs.vcf.gz", record_idx, sample_idx, 
        des = "target.chr$chr.typedOnly.maf$maf.vcf.gz", allow_multiallelic=false)

    # unphase and mask 1% entries in target file
    println("unphasing and masking entries in target file with typed snps only")
    masks = falses(p, n)
    missingprop = 0.001
    for j in 1:n, i in 1:p
        rand() < missingprop && (masks[i, j] = true)
    end
    @time mask_gt("target.chr$chr.typedOnly.maf$maf.vcf.gz", masks, 
        des="target.chr$chr.typedOnly.maf$maf.masked.vcf.gz", unphase=true)

    # compress reference file to jlso format
    max_d   = [1000]
    reffile = "ref.chr$chr.excludeTarget.vcf.gz"
    tgtfile = "target.chr$chr.typedOnly.maf$maf.masked.vcf.gz"
    H, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt = convert_ht(Bool, reffile, trans=true, save_snp_info=true, msg="importing reference data...")
    X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt = VCFTools.convert_gt(UInt8, tgtfile, trans=true, save_snp_info=true, msg = "Importing genotype file...")
    for d in max_d
        outfile = "ref.chr$chr.excludeTarget.maxd$d.jlso"
        @time compress_haplotypes(H, X, outfile, X_pos, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt, d, 0, 0.0)
    end
end
Random.seed!(2020)
data = "./chr20.1kg.phase3.v5a.vcf.gz"
maf = 0.4
chr = 20
@time filter_and_mask(data, maf, chr)


