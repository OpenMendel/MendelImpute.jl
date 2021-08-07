# MendelImpute

| **Documentation** | **Build Status** | **Code Coverage**  |
|-------------------|------------------|--------------------|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OpenMendel.github.io/MendelImpute.jl/dev/) [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenMendel.github.io/MendelImpute.jl/stable/) | [![build Actions Status](https://github.com/OpenMendel/MendelImpute.jl/workflows/CI/badge.svg)](https://github.com/OpenMendel/MendelImpute.jl/actions) [![CI (Julia nightly)](https://github.com/openmendel/mendelimpute.jl/workflows/JuliaNightly/badge.svg)](https://github.com/OpenMendel/MendelImpute.jl/actions/workflows/JuliaNightly.yml) | [![codecov](https://codecov.io/gh/OpenMendel/MendelImpute.jl/branch/master/graph/badge.svg?token=YyPqiFpIM1)](https://codecov.io/gh/OpenMendel/MendelImpute.jl) |

## Installation

Download and install [Julia](https://julialang.org/downloads/). Within Julia, copy and paste the following: 
```julia
using Pkg
pkg"add https://github.com/OpenMendel/MendelImpute.jl"
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

Our paper is on [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.10.24.353755v2). If you want to cite `MendelImpute.jl`, please cite

```
@article{mendelimpute,
    title = {{A Fast Data-Driven Method for Genotype Imputation, Phasing, and Local Ancestry Inference: MendelImpute.jl}},
    author = {Chu, Benjamin B and Sobel, Eric M and Wasiolek, Rory and Sinsheimer, Janet S and Zhou, Hua and Lange, Kenneth},
    year = {2020},
    journal={arXiv preprint DOI:10.1101/2020.10.24.353755}
}
```

## Acknowledgement

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
