using VCFTools

truth  = "target.full.vcf.gz"
mendel = "mendel.imputed.vcf.gz"
beagle = "beagle.imputed.maf0.05.vcf.gz"
mmac4  = "minimac4.result.dose.vcf.gz"
imp5   = "impute5.result.vcf.gz"

Xtrue = convert_gt(UInt8, truth, msg="importing...")
Xmendel = convert_gt(UInt8, mendel, msg="importing...")
Xbeagle = convert_gt(UInt8, beagle, msg="importing...")
Xmmac4 = convert_gt(UInt8, mmac4, msg="importing...")
Ximp5  = convert_gt(UInt8, imp5, msg="importing...")
n, p = size(Xtrue)

println("error mendel = $(sum(Xmendel .!= Xtrue) / n / p)")
println("error beagle = $(sum(Xbeagle .!= Xtrue) / n / p)")
println("error minimac4 = $(sum(Xmmac4 .!= Xtrue) / n / p)")
println("error impute5 = $(sum(Ximp5 .!= Xtrue) / n / p)")
