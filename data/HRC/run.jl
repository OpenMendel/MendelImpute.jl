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
using TimerOutputs
using LinearAlgebra

BLAS.set_num_threads(1)

Threads.nthreads()

Random.seed!(2020)
width   = 64
tgtfile = "target.chr20.typedOnly.maf0.01.masked.vcf.gz"
reffile = "ref.chr20.w$width.maf0.01.excludeTarget.jlso"
outfile = "mendel.chr20.imputed.target.vcf.gz"
@time ph = phase(tgtfile, reffile, outfile = outfile, width = width,
    dynamic_programming = false);

# import imputed result and compare with true
X_mendel = convert_gt(Float32, outfile);
X_complete = convert_gt(Float32, "target.chr20.full.vcf.gz");
n, p = size(X_mendel);
error_rate = sum(X_mendel .!= X_complete) / n / p
