using Revise
using VCFTools
using SnpArrays
using MendelImpute
using Random
using StatsBase
using CodecZlib
using ProgressMeter
using JLSO
using BenchmarkTools
using LinearAlgebra
using CSV
using DataFrames

function filter_and_mask(chr::Int)
    # filter chromosome data for unique snps
#         data = "../beagle_raw/chr$chr.1kg.phase3.v5a.vcf.gz"
#         full_record_index = .!find_duplicate_marker(data)
#         VCFTools.filter(data, full_record_index, 1:nsamples(data), 
#             des = "chr$chr.uniqueSNPs.vcf.gz", allow_multiallelic=false)

    # import VCF data with only unique SNPs
    _, vcf_sampleID, _, vcf_record_pos, _, _, _ = convert_gt(UInt8, 
        "chr$chr.uniqueSNPs.vcf.gz", save_snp_info=true, msg="importing")
    total_snps = length(vcf_record_pos)
    samples = length(vcf_sampleID)

    # generate target panel with all snps
    n = 100
    sample_idx = falses(samples)
    sample_idx[1:n] .= true
    shuffle!(sample_idx)
    VCFTools.filter("chr$chr.uniqueSNPs.vcf.gz", 1:total_snps, 
        sample_idx, des = "target.chr$chr.full.vcf.gz", allow_multiallelic=false)

    # generate reference panel without target samples
    VCFTools.filter("chr$chr.uniqueSNPs.vcf.gz", 1:total_snps, 
        .!sample_idx, des = "ref.chr$chr.excludeTarget.vcf.gz", allow_multiallelic=false)

    # generate target file with 100 samples whose snps are in the UK Biobank
    bed = CSV.read("/home/biona001/omni_chips/InfiniumOmni5-4v1-2_A1.bed", 
        skipto = 2, header=false, delim='\t', DataFrame)
    typed_snppos = bed[findall(x -> x == "chr$chr", bed[!, 1]), 3]
    match_idx = indexin(vcf_record_pos, typed_snppos)
    record_idx = falses(total_snps)
    record_idx[findall(!isnothing, match_idx)] .= true
    VCFTools.filter("chr$chr.uniqueSNPs.vcf.gz", record_idx, sample_idx, 
        des = "target.chr$chr.typedOnly.vcf.gz", allow_multiallelic=false)

    # unphase and mask 0.1% entries in target file
    p = nrecords("target.chr$chr.typedOnly.vcf.gz")
    masks = falses(p, n)
    missingprop = 0.001
    for j in 1:n, i in 1:p
        rand() < missingprop && (masks[i, j] = true)
    end
    mask_gt("target.chr$chr.typedOnly.vcf.gz", masks, 
        des="target.chr$chr.typedOnly.masked.vcf.gz", unphase=true)
end 
Random.seed!(2020)
@time filter_and_mask(10)

