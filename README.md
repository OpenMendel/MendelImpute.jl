# MendelImpute

| **Documentation** | **Build Status** | **Code Coverage**  |
|-------------------|------------------|--------------------|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://biona001.github.io/MendelImpute/dev/) | [![Build Status](https://travis-ci.com/biona001/MendelImpute.svg?branch=master)](https://travis-ci.com/github/biona001/MendelImpute) | [![Coverage Status](https://coveralls.io/repos/github/biona001/MendelImpute/badge.svg?branch=master)](https://coveralls.io/github/biona001/MendelImpute?branch=master) |

## Installation

Download and install [Julia](https://julialang.org/downloads/). Within Julia, copy and paste the following: 
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/OpenMendel/SnpArrays.jl.git"))
Pkg.add(PackageSpec(url="https://github.com/OpenMendel/VCFTools.jl.git"))
Pkg.add(PackageSpec(url="https://github.com/OpenMendel/MendelImpute.jl.git"))
```
This package supports Julia `v1.5`+.

## Documentation

+ [**Latest**](https://biona001.github.io/MendelImpute/dev/)

## Usage

Given a target genotype file (phased or unphased, may contain missing data, ending in `.vcf` or `.vcf.gz`) and a reference haplotype file (phased, no missing, ending in `.vcf`, `.vcf.gz` or `.jlso`), our software phases and imputes every SNP in the reference file to the target file. Like many other software, SNPs typed in target must all be present in the reference panel.

Also check out [VCFTools.jl](https://github.com/OpenMendel/VCFTools.jl) and [SnpArrays.jl](https://github.com/OpenMendel/SnpArrays.jl) for intuitive manipulation of VCF and PLINK files. 

## Bug Fixes and User support

If you encounter a bug or need user support, please open a new issue on Github. Please provide as much detail as possible for bug reports, ideally a sequence of reproducible code that lead to the error. 

PRs and feature requests are welcomed!

## Citation

The manuscript is still in preparation. Support us by pressing the star button on the upper right corner! 

## Acknowledgement

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
