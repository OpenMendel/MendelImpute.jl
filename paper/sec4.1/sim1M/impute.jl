using MendelImpute
using Random 
using LinearAlgebra
using VCFTools

show_error = false
warmup = true
curdir = pwd()

BLAS.set_num_threads(1)
Threads.nthreads() != 10 && error("not 10 threads!")

if warmup
    cd(normpath(MendelImpute.datadir()))
    phase("target.typedOnly.masked.vcf.gz", "ref.excludeTarget.jlso", "imputed.vcf.gz")
end

Random.seed!(2020)
cd(curdir)
tgtfile = "target.typedOnly.maf0.05.masked.vcf.gz"
reffile = "ref.excludeTarget.maxd1000.jlso"
outfile = "mendel.imputed.jlso"
phase(tgtfile, reffile, outfile)

if show_error
    X_complete = convert_gt(UInt8, "target.full.vcf.gz")
    n, p = size(X_complete)
    X_mendel = convert_gt(UInt8, outfile)
    println("error overall = $(sum(X_mendel .!= X_complete) / n / p) \n")
end
