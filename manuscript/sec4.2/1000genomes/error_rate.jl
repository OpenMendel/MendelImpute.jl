using VCFTools

chr = 20
truth  = "target.chr$chr.full.vcf.gz"
mendel = "mendel.imputed.chr$chr.vcf.gz"
beagle = "beagle.imputed.chr$chr.vcf.gz"
mmac4  = "minimac4.chr$chr.result.dose.vcf.gz"

isfile(truth) && isfile(mendel) && isfile(beagle) && isfile(mmac4) || error("file doens't exist!")

Xtrue = convert_gt(UInt8, truth, msg="importing Xtrue")
Xmendel = convert_gt(UInt8, mendel, msg="importing Xmendel")
Xbeagle = convert_gt(UInt8, beagle, msg="importing Xbeagle")
Xmmac4 = convert_gt(UInt8, mmac4, msg="importing Xmmac")
n, p = size(Xtrue)

# need to import each chunk of impute5 result separately
chunks = chr == 20 ? 3 : 7
Ximp5 = Vector{Matrix{UInt8}}(undef, chunks)
for i in 1:chunks
    Ximp5[i] = convert_gt(UInt8, "impute5.chr$chr.result.chunk$i.vcf.gz", msg="importing Ximp5 $i / $chunks")
end
Ximp5 = hcat(Ximp5...)

println("error mendel = $(sum(Xmendel .!= Xtrue) / n / p)")
println("error beagle = $(sum(Xbeagle .!= Xtrue) / n / p)")
println("error minimac4 = $(sum(Xmmac4 .!= Xtrue) / n / p)")
println("error impute5 = $(sum(Ximp5 .!= Xtrue) / n / p)")
