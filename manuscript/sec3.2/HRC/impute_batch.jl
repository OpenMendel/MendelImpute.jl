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

BLAS.set_num_threads(1)
Threads.nthreads() != 10 && error("not 10 threads!")

Random.seed!(2020)
function run()
    chr = 20
    maf = 0.01
    max_d = 1000
    X_complete = convert_gt(UInt8, "target.chr$chr.full.vcf.gz")
    n, p = size(X_complete)
    for overlap in [0.25, 0.5, 1.0]
	tgtfile = "target.chr$chr.typedOnly.masked.vcf.gz"
	reffile = "ref.chr$chr.maxd$max_d.overlap$overlap.excludeTarget.jlso"
	outfile = "mendel.imputed.chr$chr.vcf.gz"

	ph = phase(tgtfile, reffile, outfile=outfile, impute=true, max_d=max_d, phase=true)
        X_mendel = convert_gt(UInt8, outfile)
        println("overlap $overlap error = $(sum(X_mendel .!= X_complete) / n / p) \n")
        flush(stdout)
    end
end
run()
