web = "http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/"
vcf = "1kg.phase3.v5a.vcf.gz"
tbi = "1kg.phase3.v5a.vcf.gz.tbi"

for chr in 1:22
    println("downloading chromosome $chr data")
    vcffile = web * "chr$chr." * vcf
    tbifile = web * "chr$chr." * tbi
    isfile(vcffile) || download(vcffile, "chr$chr." * vcf)
    isfile(tbifile) || download(tbifile, "chr$chr." * tbi)
end

vcffile = web * "chrX." * vcf
tbifile = web * "chrX." * tbi
isfile(vcffile) || download(vcffile, "chrX." * vcf)
isfile(tbifile) || download(tbifile, "chrX." * tbi)
