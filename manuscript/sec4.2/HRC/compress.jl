using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using StatsBase
using CodecZlib
using ProgressMeter
using JLD2, FileIO, JLSO
using BenchmarkTools
using GroupSlices
using LinearAlgebra
using ThreadPools

function compress()
        max_d = [1000]
	chr = 20
	#minwidths = [50, 100, 250, 500]
	#overlaps = [0.0]
	minwidths = [0]
	overlap = 0.0

        tgtfile = "target.chr$chr.typedOnly.masked.vcf.gz"
        reffile = "ref.chr$chr.excludeTarget.vcf.gz"
	X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt = VCFTools.convert_gt(UInt8, tgtfile, trans=true, save_snp_info=true, msg = "Importing genotype file...")
        H, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt = convert_ht(Bool, reffile, trans=true, save_snp_info=true, msg="importing reference data...")
        for d in max_d, minwidth in minwidths
            outfile = "ref.chr$chr.maxd$d.excludeTarget.jlso"
            @time compress_haplotypes(H, X, outfile, X_pos, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt, d, minwidth, overlap)
        end
end
compress()
