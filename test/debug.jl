using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using SparseArrays
using JLD2, FileIO, JLSO
using ProgressMeter
using GroupSlices
using ThreadPools

Random.seed!(2020)
width   = 512
tgtfile = "./compare2/target.typedOnly.maf0.01.masked.vcf.gz"
reffile = "./compare2/ref.excludeTarget.w$width.jlso"
outfile = "./compare2/mendel.imputed.vcf.gz"
@time ph = phase(tgtfile, reffile, outfile = outfile, width = width,
    dynamic_programming = false);
