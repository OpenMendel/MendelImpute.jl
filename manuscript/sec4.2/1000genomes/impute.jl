using Revise
using VCFTools
using MendelImpute
using Random
using LinearAlgebra

chr = 10
d = 1000
show_error = false
warmup = true
curdir = pwd()

BLAS.set_num_threads(1)
Threads.nthreads() != 10 && error("not 10 threads!")

if warmup
    cd(normpath(MendelImpute.datadir()))
    #compress_haplotypes("ref.excludeTarget.vcf.gz", "target.typedOnly.masked.vcf.gz", "ref.excludeTarget.jlso")
    phase("target.typedOnly.masked.vcf.gz", "ref.excludeTarget.jlso", "imputed.vcf.gz")
end

Random.seed!(2020)
cd(curdir)
tgtfile = "target.chr$chr.typedOnly.masked.vcf.gz"
reffile = "ref.chr$chr.excludeTarget.jlso"
outfile = "mendel.imputed.chr$chr.jlso"
ph = phase(tgtfile, reffile, outfile)

if show_error
    X_complete = convert_gt(UInt8, "target.chr$chr.full.vcf.gz")
    n, p = size(X_complete)
    X_mendel = convert_gt(UInt8, outfile)
    println("error overall = $(sum(X_mendel .!= X_complete) / n / p) \n")
end
