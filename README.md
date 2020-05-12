# MendelImpute

[![Build Status](https://travis-ci.com/biona001/MendelImpute.svg?branch=master)](https://travis-ci.com/github/biona001/MendelImpute) [![Coverage Status](https://coveralls.io/repos/github/biona001/MendelImpute/badge.svg?branch=master)](https://coveralls.io/github/biona001/MendelImpute?branch=master)

## Installation

Start Julia, press ] to enter package manager mode, and type:
```julia
(v1.3) pkg> add https://github.com/OpenMendel/VCFTools.jl
(v1.3) pkg> add https://github.com/biona001/MendelImpute
```

## Usage

Our software takes a data-mining approach for genotype imputation, contrary to HMM or low rank approximation methods. Given a target genotype file (phased or unphased and may contain missing data) and a reference haplotype file (phased, no missing), our software phases and imputes every SNP in the reference file to the target file. Observed data remains unchanged.

Manuscript coming soon!

## Examples

Please impute one chromosome at a time. 

```julia
tgtfile = "target.vcf.gz"           # Target imputation file name. Doesn't have to be phased. 
reffile = "ref.vcf.gz"              # Reference haplotype file name. Genotypes must be phased. 
outfile = "mendel.imputed.vcf.gz"   # Output file name. All output genotypes will be phased. 
width = 400                         # Number of SNPs per imputation window. (default 400)
@time phase(tgtfile, reffile, outfile = outfile, impute=true, width = width)
```
Example output:
```julia
Importing genotype file...100%|█████████████████████████| Time: 0:00:05
Importing reference haplotype files...100%|█████████████| Time: 0:03:07
Computing optimal haplotype pairs...100%|███████████████| Time: 0:03:23
Merging breakpoints...100%|█████████████████████████████| Time: 0:01:14
Writing to file...100%|█████████████████████████████████| Time: 0:00:17
517.542827 seconds (3.07 G allocations: 288.416 GiB, 7.70% gc time)
```

Check out our [Jupyter notebook examples](https://github.com/biona001/MendelImpute/tree/master/data/1000_genome_phase3_v5/filtered) (click on any `.pynb` link)! The result above is from `chrom22` notebook ran on 8 threads, with 644939 reference SNPs, 100000 target SNPs, 250 samples, and 4508 haplotypes.

## Options

These options are usable with the `phase` function.

- `outfile`: output filename. Output genotypes will be phased with no missing data.
- `impute`: If `true`, untyped SNPs will be imputed, otherwise only missing snps in `tgtfile` will be imputed. (default `false`)
- `width`: number of SNPs (markers) in each sliding window. (default `400`)
- `flankwidth`: Number of SNPs flanking the sliding window (defaults to 10% of `width`)
- `fast_method`: If `true`, will use window-by-window intersection for phasing. If `false`, phasing uses dynamic progrmaming.  (default `false`)
- `unique_only`: If `true`, will phase and impute using only unique haplotypes in each window. Usually faster than `fast_method=true`.  (default `false`)

**Note for multithreading**: The default number of threads is 1. To change this, type `export JULIA_NUM_THREADS=4` in your terminal *before starting julia*. We recommend setting number of threads equal to number of physical CPU cores. 

## Package features

- Built-in support for `.vcf`, `.vcf.gz` and PLINK (coming soon) files.
- Out-of-the-box multi-threaded parallelism (before starting julia, type `export JULIA_NUM_THREADS=4`)
- Impute dosage data (genotype is any real number in [0, 2]) using a haplotype reference panel (coming soon)
- Intuitive manipulation of genotype files via `VCFTools.jl` and `SnpArrays.jl`
- Some simple simulation routines for generating haplotype and phase/unphased target files. 

## Bug Fixes and User support

If you encounter a bug or need user support, please open a new issue on Github. Please provide as much detail as possible for bug reports, ideally a sequence of reproducible code that lead to the error. 

PRs and feature requests are welcomed!

## Citation

The manuscript is still in preparation. Support us by pressing the star button on the upper right corner! 
