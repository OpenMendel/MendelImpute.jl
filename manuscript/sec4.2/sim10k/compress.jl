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
	minwidth = 0
	overlap = 0.0

        tgtfile = "target.typedOnly.maf0.05.masked.vcf.gz"
        reffile = "ref.excludeTarget.vcf.gz"
        outfile = "ref.excludeTarget.maxd1000.jlso"
        @time compress_haplotypes(reffile, tgtfile, outfile, max_d, minwidth, overlap)
end
compress()
