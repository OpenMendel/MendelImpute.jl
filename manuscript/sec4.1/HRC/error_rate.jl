using VCFTools

chr = 20
truth  = "target.chr$chr.full.vcf.gz"
mendel = "mendel.imputed.chr$chr.vcf.gz"
beagle = "beagle.chr$chr.imputed.vcf.gz"
mmac4  = "minimac4.chr$chr.result.dose.vcf.gz"

isfile(truth) && isfile(mendel) && isfile(beagle) && isfile(mmac4) || error("file doens't exist!")

Xtrue = convert_gt(UInt8, truth)
Xmendel = convert_gt(UInt8, mendel)
Xbeagle = convert_gt(UInt8, beagle)
Xmmac4 = convert_gt(UInt8, mmac4)
n, p = size(Xtrue)

println("error mendel = $(sum(Xmendel .!= Xtrue) / n / p)")
println("error beagle = $(sum(Xbeagle .!= Xtrue) / n / p)")
println("error minimac4 = $(sum(Xmmac4 .!= Xtrue) / n / p)")
