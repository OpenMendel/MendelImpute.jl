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
using Profile

chr = 20
maf = 0.01
width = 64

BLAS.set_num_threads(1)
Threads.nthreads() != 10 && error("not 10 threads!")
Profile.init(n = 10^7, delay = 1.0)

Random.seed!(2020)
tgtfile = "target.chr$chr.typedOnly.maf$maf.masked.vcf.gz"
reffile = "./single_precision/ref.chr$chr.w$width.maf$maf.excludeTarget.jlso"
outfile = "mendel.imputed.dp$width.maf$maf.vcf.gz"

phase(tgtfile, reffile, outfile=outfile, impute=true, width=width, 
    dynamic_programming = false)
Profile.clear_malloc_data()
phase(tgtfile, reffile, outfile=outfile, impute=true, width=width,
    dynamic_programming = false)

rm(outfile, force=true)
