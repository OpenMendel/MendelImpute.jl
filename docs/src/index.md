# MendelImpute.jl 

*Fast genotype imputation, phasing, and admixture estimation!*

`MendelImpute.jl` is the fastest and least memory-consuming software for phasing and genotype imputation, as of 2020. It is also capable of ancestry estimation.

## Package Features

- Built-in support for imputing genotypes stored in [VCF files](https://samtools.github.io/hts-specs/VCFv4.3.pdf) (`.vcf`, `.vcf.gz`) or [PLINK files](https://www.cog-genomics.org/plink2/formats#bed).
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
    "man/Phasing+and+Imputation.md"
    "man/performance.md"
    "man/painting.md"
    "man/api.md"
]
Depth = 2
```