using Revise
using VCFTools
using MendelImpute
using Random
using LinearAlgebra

chr = 20
d = 1000
show_error = false
warmup = true
curdir = pwd()
maf = chr == 10 ? 0.25 : 0.01

BLAS.set_num_threads(1)
Threads.nthreads() != 10 && error("not 10 threads!")

if warmup
    cd(normpath(MendelImpute.datadir()))
    phase("target.typedOnly.masked.vcf.gz", "ref.excludeTarget.jlso", "imputed.vcf.gz")
end

Random.seed!(2020)
cd(curdir)
tgtfile = "target.chr$chr.typedOnly.maf$maf.masked.vcf.gz"
reffile = "ref.chr$chr.excludeTarget.maxd$d.jlso"
outfile = "mendel.imputed.chr$chr.vcf.gz"
ph = phase(tgtfile, reffile, outfile)

if show_error
    X_complete = convert_gt(UInt8, "target.chr$chr.full.vcf.gz")
    n, p = size(X_complete)
    X_mendel = convert_gt(UInt8, outfile)
    println("error overall = $(sum(X_mendel .!= X_complete) / n / p) \n")
end
