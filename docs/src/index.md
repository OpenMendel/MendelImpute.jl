# MendelImpute.jl 

*Fast genotype imputation, phasing, and admixture estimation!*

`MendelImpute.jl` is the fastest and most memory-efficient software for phasing and genotype imputation, as of 2021. It is also capable of estimating local and global admixture coefficients using a reference haplotype panel.

Given a target genotype file (phased or unphased and may contain missing data) and a reference haplotype file (phased, no missing), our software imputes every SNP in the reference file to the target file, outputing phased or unphased genotypes. Like many other software, SNPs typed in target must all be present in the reference panel. 

## Package Features

- Phasing and imputation with respect to a reference haplotype panel
- Imputation on dosage data, phasing without imputation, imputation without phasing
- Built-in support for [VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf) (`.vcf`, `.vcf.gz`), [PLINK](https://www.cog-genomics.org/plink2/formats#bed), and [BGEN](https://www.well.ox.ac.uk/~gav/bgen_format/)  (`.bgen`, currently experimental) files
- Out-of-the-box multithreaded (shared memory) parallelism. 
- Admixture estimation, with code examples to make pretty plots!
- Ultra-compressed file for phased genotypes.

## Installation

Download and install [Julia](https://julialang.org/downloads/). Within Julia, copy and paste the following: 
```julia
using Pkg
pkg"add https://github.com/OpenMendel/SnpArrays.jl"
pkg"add https://github.com/OpenMendel/VCFTools.jl"
pkg"add https://github.com/OpenMendel/BGEN.jl"
pkg"add https://github.com/OpenMendel/MendelImpute.jl"
```
This package supports Julia `v1.6`+.

## Manual Outline

```@contents
Pages = [
    "man/Phasing_and_Imputation.md"
    "man/performance.md"
    "man/painting.md"
    "man/ultra+compress.md"
    "man/script.md"
    "man/api.md"
]
Depth = 2
```
