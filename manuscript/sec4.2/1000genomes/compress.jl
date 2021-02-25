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
        max_d = 1000
	chr = 10
	minwidth = 0
	overlap = 0.0

        tgtfile = "target.chr$chr.typedOnly.masked.vcf.gz"
        reffile = "ref.chr$chr.excludeTarget.vcf.gz"
        outfile = "ref.chr$chr.excludeTarget.jlso"
        @time compress_haplotypes(reffile, tgtfile, outfile, max_d, minwidth, overlap)
end
compress()
