# MendelImpute

| **Documentation** | **Build Status** | **Code Coverage**  |
|-------------------|------------------|--------------------|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OpenMendel.github.io/MendelImpute.jl/dev/) | [![Build Status](https://travis-ci.com/OpenMendel/MendelImpute.jl.svg?branch=master)](https://travis-ci.com/github/OpenMendel/MendelImpute.jl) | [![Coverage Status](https://coveralls.io/repos/github/OpenMendel/MendelImpute.jl/badge.svg?branch=master)](https://coveralls.io/github/OpenMendel/MendelImpute.jl?branch=master) |

## Installation

Download and install [Julia](https://julialang.org/downloads/). Within Julia, copy and paste the following: 
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/OpenMendel/SnpArrays.jl.git"))
Pkg.add(PackageSpec(url="https://github.com/OpenMendel/VCFTools.jl.git"))
Pkg.add(PackageSpec(url="https://github.com/OpenMendel/MendelImpute.jl.git"))
```
This package supports Julia `v1.5`+.

Example run:

```julia
# first compress reference haplotypes to .jlso forma
using MendelImpute                         # load package
cd(normpath(MendelImpute.datadir())        # change to data directory
reffile = "ref.excludeTarget.vcf.gz"       # specify reference VCF file
tgtfile = "target.typedOnly.masked.vcf.gz" # specify target VCF file (GWAS file)
outfile = "ref.excludeTarget.jlso"         # output file name (end in .jlso)
compress_haplotypes(reffile, tgtfile, outfile)

# phase & impute
tgtfile = "target.typedOnly.masked.vcf.gz" # target VCF file (GWAS file)
reffile = "ref.excludeTarget.jlso"         # compressed ref file
outfile = "imputed.vcf.gz"                 # output file name (phased & imputed GWAS data)
phase(tgtfile, reffile, outfile = outfile)
```

## Documentation

+ [**Latest**](https://OpenMendel.github.io/MendelImpute.jl/dev/)

## Bug Fixes and User support

If you encounter a bug or need user support, please open a new issue on Github. Please provide as much detail as possible for bug reports, ideally a sequence of reproducible code that lead to the error. 

PRs and feature requests are welcomed!

## Citation

The manuscript is still in preparation. Support us by pressing the star button on the upper right corner! 

## Acknowledgement

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
