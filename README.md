# MendelImpute

| **Documentation** | **Build Status** | **Code Coverage**  |
|-------------------|------------------|--------------------|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OpenMendel.github.io/MendelImpute.jl/dev/) [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenMendel.github.io/MendelImpute.jl/stable/) | [![build Actions Status](https://github.com/OpenMendel/MendelImpute.jl/workflows/CI/badge.svg)](https://github.com/OpenMendel/MendelImpute.jl/actions) [![CI (Julia nightly)](https://github.com/openmendel/mendelimpute.jl/workflows/JuliaNightly/badge.svg)](https://github.com/OpenMendel/MendelImpute.jl/actions/workflows/JuliaNightly.yml) | [![codecov](https://codecov.io/gh/OpenMendel/MendelImpute.jl/branch/master/graph/badge.svg?token=YyPqiFpIM1)](https://codecov.io/gh/OpenMendel/MendelImpute.jl) |

## Installation

Download and install [Julia](https://julialang.org/downloads/). Within Julia, copy and paste the following: 
```julia
using Pkg
pkg"add MendelImpute"
```
This package supports Julia `v1.6`+. 

## Documentation

+ [**Latest**](https://OpenMendel.github.io/MendelImpute.jl/dev/)
+ [**Stable**](https://OpenMendel.github.io/MendelImpute.jl/stable/)

## Example run:

The following uses data under the `data/` directory.

```julia
# load package & cd to data directory
using MendelImpute
cd(normpath(MendelImpute.datadir()))

# compress reference haplotypes from .vcf.gz to .jlso format
reffile = "ref.excludeTarget.vcf.gz"       # reference VCF file
tgtfile = "target.typedOnly.masked.vcf.gz" # target VCF file (GWAS file)
outfile = "ref.excludeTarget.jlso"         # output file name (end in .jlso)
@time compress_haplotypes(reffile, tgtfile, outfile)

# phase & impute (note: 2nd run will be much faster because code is compiled)
tgtfile = "target.typedOnly.masked.vcf.gz" # target VCF file (GWAS file)
reffile = "ref.excludeTarget.jlso"         # compressed reference file
outfile = "imputed.vcf.gz"                 # output file name
@time phase(tgtfile, reffile, outfile);

# check error rate (since data was simulated)
using VCFTools
Ximputed = convert_gt(Float64, "imputed.vcf.gz")  # imputed genotypes
Xtrue = convert_gt(Float64, "target.full.vcf.gz") # true genotypes
m, n = size(Xtrue) # matrix dimensions
error_rate = sum(Xtrue .!= Ximputed) / m / n
```

We also support PLINK binary files (`.bed/.bim/.fam`) via [SnpArrays.jl](https://github.com/OpenMendel/SnpArrays.jl) and BGEN files `.bgen` via [BGEN.jl](https://github.com/OpenMendel/BGEN.jl). 

For more realistic example, see [detailed example in documentation](https://openmendel.github.io/MendelImpute.jl/dev/man/Phasing_and_Imputation/#Detailed-Example)

## Bug Fixes and User support

If you encounter a bug or need user support, please open a new issue on Github. Please provide as much detail as possible for bug reports, ideally a sequence of reproducible code that lead to the error. 

PRs and feature requests are welcomed!

## Citation

If you use this analysis package in your research, please cite the following references in resulting publications:

*Chu BB, Sobel EM, Wasiolek R, Ko S, Sinsheimer JS, Zhou H, Lange K. A fast Data-Driven method for genotype imputation, phasing, and local ancestry inference: MendelImpute.jl. Bioinformatics. 2021 Jul 21;37(24):4756–63. doi: 10.1093/bioinformatics/btab489. Epub ahead of print. PMID: 34289008; PMCID: [PMC8665755](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8665755/).*

*Zhou H, Sinsheimer JS, Bates DM, Chu BB, German CA, Ji SS, Keys KL, Kim J, Ko S, Mosher GD, Papp JC, Sobel EM, Zhai J, Zhou JJ, Lange K. OPENMENDEL: a cooperative programming project for statistical genetics. Hum Genet. 2020 Jan;139(1):61-71. doi: 10.1007/s00439-019-02001-z. Epub 2019 Mar 26. PMID: 30915546; PMCID: [PMC6763373](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6763373/).*

## Acknowledgement

This project has been supported by the National Institutes of Health under awards R01GM053275, R01HG006139, R25GM103774, and 1R25HG011845.
