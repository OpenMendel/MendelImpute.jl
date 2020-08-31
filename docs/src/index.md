# MendelImpute.jl 

*Fast genotype imputation, phasing, and admixture estimation!*

`MendelImpute.jl` is the fastest and least memory-consuming software for phasing and genotype imputation, as of 2020. It is also capable of ancestry estimation.

## Package Features

- Built-in support for imputing VCF (`.vcf`, `.vcf.gz`) and PLINK files.
- Out-of-the-box multithreaded (shared memory) parallelism. 
- Admixture estimation, with code examples to make pretty plots!
- Ultra-compressed file for phased genotypes.
- Imputation on dosage data.

Given a target genotype file (phased or unphased and may contain missing data) and a reference haplotype file (phased, no missing), our software imputes every SNP in the reference file to the target file, outputing phased or unphased genotypes. Like many other software, SNPs typed in target must all be present in the reference panel. 

## Installation

Download and install [Julia](https://julialang.org/downloads/). Within Julia, copy and paste the following: 
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/OpenMendel/SnpArrays.jl.git"))
Pkg.add(PackageSpec(url="https://github.com/OpenMendel/VCFTools.jl.git"))
Pkg.add(PackageSpec(url="https://github.com/OpenMendel/MendelImpute.jl.git"))
```
This package supports Julia `v1.5`+.

## Manual Outline

```@contents
Pages = [
    "man/api.md"
]
Depth = 2
```