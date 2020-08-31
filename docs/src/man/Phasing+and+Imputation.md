
# Preparing Target Data

MendelImpute accepts VCF, compressed VCF, and PLINK files. Please make sure the following are true:

+ VCF file ends in `.vcf` or `.vcf.gz`
+ For PLINK files, all trios (`.bim`, `.bed`, `.fam`) are present in the same directory
+ Each file contains only 1 chromosome
+ Every record (SNP) is present in the reference panel. If this is untrue, you must [match markers in 2 VCF files](https://openmendel.github.io/VCFTools.jl/dev/man/conformgt/). 
+ The position of every SNP is unique (so multiallelic markers should be excluded instead of split)

If the last criteria is not met, our code may or may not work. File an issue to let us know.

# Preparing Reference Haplotype Panel

Reference panels must be compressed into `.jlso` format first using the [compress_haplotypes](https://biona001.github.io/MendelImpute/dev/man/api/#MendelImpute.compress_haplotypes) function. One must specify `d`: the maximum number of unique haplotypes per window. Larger `d` slows down computation, but increases accuracy. For most purposes, we recommend $d \approx 1000$. A larger `d` may be needed for TOPMed or HRC data. 

# Detailed Example

We use the [1000 genomes chromosome 22](http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/) as an example. 

## Step 1: generating realistic reference and target data 

First we generate a reference panel and imputation target based on the 1000 genomes data. More specifically, 
+ The first 100 samples are used as imputation targets, where
    - 100k SNPs with minor allele frequency $\ge 0.05$ are randomly selected to be the typed positions. 
    - 0.1% of typed SNPs are masked (mimicking GWAS errors)
    - Genotypes are unphased
+ The remaining 2404 samples are used as reference haplotypes. 
+ SNPs with duplicate positions are filtered out.
+ All multiallelic markers are filtered out.

**Instruction: execute the code below in a Julia session or a Jupyter notebook:**


```julia
# load necessary packages in Julia
using MendelImpute
using VCFTools
using Random

# set random seed for reproducibility
Random.seed!(2020)

# download example data 
data = "chr22.1kg.phase3.v5a.vcf.gz"
if !isfile(data) 
    download("http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr22.1kg.phase3.v5a.vcf.gz")
end

# remove SNPs with the same positions, keep all samples, save result into new file
SNPs_to_keep = .!find_duplicate_marker(data) 
@time VCFTools.filter(data, SNPs_to_keep, 1:nsamples(data), des = "chr22.uniqueSNPs.vcf.gz")

# summarize data
total_snps, samples, _, _, _, maf_by_record, _ = gtstats("chr22.uniqueSNPs.vcf.gz")

# generate target file with 100 samples and 100k snps with maf>0.05
n = 100
p = 100000
record_idx = falses(total_snps)
large_maf = findall(x -> x > 0.05, maf_by_record)  
Random.shuffle!(large_maf)
record_idx[large_maf[1:p]] .= true
sample_idx = falses(samples)
sample_idx[1:n] .= true
Random.shuffle!(sample_idx)
@time VCFTools.filter("chr22.uniqueSNPs.vcf.gz", record_idx, sample_idx, 
    des = "target.chr22.typedOnly.vcf.gz", allow_multiallelic=false)

# unphase and mask 0.1% entries in target file
masks = falses(p, n)
missingprop = 0.001
for j in 1:n, i in 1:p
    rand() < missingprop && (masks[i, j] = true)
end
@time mask_gt("target.chr22.typedOnly.vcf.gz", masks, 
    des="target.chr22.typedOnly.masked.vcf.gz", unphase=true)

# generate target panel with all snps (this file contains true phase and genotypes)
@time VCFTools.filter("chr22.uniqueSNPs.vcf.gz", 1:total_snps, 
    sample_idx, des = "target.chr22.full.vcf.gz", allow_multiallelic=false)

# generate reference panel with 2404 samples
@time VCFTools.filter("chr22.uniqueSNPs.vcf.gz", 1:total_snps, .!sample_idx, 
    des = "ref.chr22.excludeTarget.vcf.gz", allow_multiallelic=false)
```

    ┌ Info: Precompiling MendelImpute [e47305d1-6a61-5370-bc5d-77554d143183]
    └ @ Base loading.jl:1278
    [32mfinding duplicate markers...100%|███████████████████████| Time: 0:03:56[39m
    [32mfiltering vcf file...100%|██████████████████████████████| Time: 0:04:46[39m


    292.131527 seconds (3.20 G allocations: 301.789 GiB, 7.89% gc time)


    [32mProgress: 100%|█████████████████████████████████████████| Time: 0:04:02[39m
    [32mfiltering vcf file...100%|██████████████████████████████| Time: 0:03:59[39m


    244.425505 seconds (3.18 G allocations: 301.694 GiB, 9.69% gc time)
      1.935526 seconds (20.00 M allocations: 1.491 GiB, 6.33% gc time)


    [32mfiltering vcf file...100%|██████████████████████████████| Time: 0:04:10[39m


    255.505399 seconds (3.27 G allocations: 317.749 GiB, 9.95% gc time)


    [32mfiltering vcf file...100%|██████████████████████████████| Time: 0:07:27[39m


    453.383147 seconds (6.16 G allocations: 566.535 GiB, 10.16% gc time)


### Output explanation:

You just generated reference and target VCF files:

+ `ref.chr22.excludeTarget.vcf.gz`: Reference haplotype panel with 2404 samples
+ `target.chr22.typedOnly.masked.vcf.gz`: Imputation target file containing 100 samples at 100k SNPs. All genotypes are unphased and contains 0.1% missing data. 

You also generated/downloaded:

+ `chr22.1kg.phase3.v5a.vcf.gz`: The original chromosome 22 data downloaded from Beagle's website.
+ `chr22.uniqueSNPs.vcf.gz`: This is the original chromosome 22 data excluding duplicate records (SNPs) by checking marker positions. The first SNP is included but all subsequent SNPs are removed. 
+ `target.chr22.full.vcf.gz`: The complete data for imputation target, used for checking imputation accuracy. All genotypes are phased and non-missing. 
+ `target.chr22.typedOnly.vcf.gz`: Complete target data on just the typed SNPs. All genotypes are phased and non-missing. Just by-producted for generating other files; not used for anything downstream.

## Step 2: generating `.jlso` compressed reference panel

MendelImpute requires one to pre-process the reference panel for faster reading. This is achieved via the [compress_haplotypes](https://biona001.github.io/MendelImpute/dev/man/api/#MendelImpute.compress_haplotypes) function.


```julia
max_d = 1000 # maximum number of unique haplotypes per window
reffile = "ref.chr22.excludeTarget.vcf.gz"
tgtfile = "target.chr22.typedOnly.masked.vcf.gz"
outfile = "ref.chr22.maxd1000.excludeTarget.jlso"
@time compress_haplotypes(reffile, tgtfile, outfile, max_d)
```

    [32mimporting reference data...100%|████████████████████████| Time: 0:01:55[39m


    284.456419 seconds (2.08 G allocations: 209.110 GiB, 10.89% gc time)


## Step 3: Run imputation and phasing

The main function is the [phase](https://biona001.github.io/MendelImpute/dev/man/api/#MendelImpute.phase) function. The code below runs it in a single thread:


```julia
reffile = "ref.chr22.maxd1000.excludeTarget.jlso" # jlso reference file
tgtfile = "target.chr22.typedOnly.masked.vcf.gz"  # target genotype file
outfile = "mendel.imputed.chr22.vcf.gz"           # output file name
d       = 1000 # this should be the equal to max_d in previous code block
phase(tgtfile, reffile; outfile=outfile, max_d = d);
```

    Importing reference haplotype data...


    [32mComputing optimal haplotypes...100%|████████████████████| Time: 0:00:21[39m


    Total windows = 1634, averaging ~ 508 unique haplotypes per window.
    
    Timings: 
        Data import                     = 11.719 seconds
            import target data             = 2.23359 seconds
            import compressed haplotypes   = 9.48538 seconds
        Computing haplotype pair        = 22.1275 seconds
            BLAS3 mul! to get M and N      = 0.978357 seconds per thread
            haplopair search               = 17.6055 seconds per thread
            initializing missing           = 0.092384 seconds per thread
            allocating and viewing         = 0.300453 seconds per thread
            index conversion               = 0.00979729 seconds per thread
        Phasing by win-win intersection = 1.45252 seconds
            Window-by-window intersection  = 0.571672 seconds per thread
            Breakpoint search              = 0.319126 seconds per thread
            Recording result               = 0.0618488 seconds per thread
        Imputation                     = 3.24156 seconds
            Imputing missing               = 0.149169 seconds
            Writing to file                = 3.09239 seconds
    
        Total time                      = 38.6959 seconds
    


Inputs after the first `;` are all optional, which is why an equal sign is needed. If left not specified, MendelImpute uses the default values. A list of optional inputs can be found in the [phase( ) API](https://biona001.github.io/MendelImpute/dev/man/api/#MendelImpute.phase). The second `;` hides the output, or else the screen will be too jammed. 

!!! note

    To run MendelImpute in parallel, type `export JULIA_NUM_THREADS=4` **before** starting Julia. For optimal performance, set threads equal to the number of physical cores on your machine. Verify the Julia session is running is parallel by executing `Threads.nthreads()` in Julia. **Finally**, set the number of BLAS threads to be 1 by `using LinearAlgebra; BLAS.set_num_threads(1)`, to avoid oversubscription. 

## Step 4: (only for simulated data) check imputation accuracy

Since we simulated data, we can check imputation accuracy.


```julia
X_truth  = convert_gt(Float64, "target.chr22.full.vcf.gz")    # import true genotypes
X_mendel = convert_gt(Float64, "mendel.imputed.chr22.vcf.gz") # import imputed genotypes
n, p = size(X_mendel)
println("error overall = $(sum(X_mendel .!= X_truth) / n / p)")
```

    error overall = 0.005397602533930585


# Run MendelImpute as script

If you don't want to run `MendelImpute.jl` in a Julia session (e.g. you want to run batch jobs on a cluster), you can do so by putting the code above in a Julia file. For instance, create a file called `impute.jl` which contains:

```julia
# place these code in a file called impute.jl
using MendelImpute, VCFTools
reffile = ARGS[1] # first command line argument
tgtfile = ARGS[2] # second command line argument
phase(tgtfile, reffile; outfile="mendel.imputed.chr22.vcf.gz", max_d = 1000)
```

Then in the terminal/command-line, you can do
```
julia impute.jl ref.chr22.maxd1000.excludeTarget.jlso target.chr22.typedOnly.masked.vcf.gz
```

