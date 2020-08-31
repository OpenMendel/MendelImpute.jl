
# Ultra-compressed format

One can optionally save/load ultra-compressed phased genotypes after imputation. Ultra-compression is nothing fancy. Instead of converting haplotype segments into genotypes, this protocol simply saves the starting position and the correct haplotype label. We put this result into our own data structure, and saving/loading is achieved by the [JLSO package](https://github.com/invenia/JLSO.jl). 

## Saving

Appending `.jlso` to the output file name will signal MendelImpute to save data in ultra-compressed format. For admixture estimation, we strongly recommend one to save in `.jlso` format.


```julia
# first load all necessary packages
using MendelImpute
using VCFTools

# compute each person's phase information
tgtfile = "target.chr22.typedOnly.masked.vcf.gz"
reffile = "ref.chr22.maxd1000.excludeTarget.jlso"
outfile = "mendel.imputed.jlso" # output file name ends in jlso!
@time phaseinfo = phase(tgtfile, reffile, outfile=outfile, impute=true, max_d=1000);
```

    Number of threads = 1
    Importing reference haplotype data...


    [32mComputing optimal haplotypes...100%|████████████████████| Time: 0:00:22[39m


    Total windows = 1634, averaging ~ 508 unique haplotypes per window.
    
    Timings: 
        Data import                     = 14.2706 seconds
            import target data             = 3.27572 seconds
            import compressed haplotypes   = 10.9948 seconds
        Computing haplotype pair        = 22.7159 seconds
            BLAS3 mul! to get M and N      = 1.00284 seconds per thread
            haplopair search               = 18.0964 seconds per thread
            initializing missing           = 0.095897 seconds per thread
            allocating and viewing         = 0.310543 seconds per thread
            index conversion               = 0.00992562 seconds per thread
        Phasing by win-win intersection = 1.30819 seconds
            Window-by-window intersection  = 0.525698 seconds per thread
            Breakpoint search              = 0.275337 seconds per thread
            Recording result               = 0.0496237 seconds per thread
        Imputation                     = 2.84315 seconds
            Imputing missing               = 0.0505533 seconds
            Writing to file                = 2.7926 seconds
    
        Total time                      = 41.2741 seconds
    
     58.892784 seconds (137.79 M allocations: 7.481 GiB, 5.92% gc time)


The object saved to `mendel.imputed.jlso` is literally the `phaseinfo` variable. We can inspect its element:


```julia
# look at sample 1's haplotype segments
haplotype_labels = phaseinfo[1].strand1.haplotypelabel # strand1
haplotype_start = phaseinfo[1].strand1.start # strand1
[haplotype_start haplotype_labels]
```




    547×2 Array{Int64,2}:
          1  4119
        236   887
        423   272
        622    12
        754   124
        802     4
        824    24
        968  1282
       1125  1741
       1202  4543
       1691  1198
       3031    22
       3521    18
          ⋮  
     411702   877
     412362    74
     413734  3849
     413868   248
     414456    31
     414552  3187
     414989  4481
     415807     5
     416108   143
     416353  1276
     416844    71
     417084   311



## Loading

The function [convert_compressed](https://biona001.github.io/MendelImpute/dev/man/api/#MendelImpute.convert_compressed) will load the ultra-compressed data into genotype matrices and the original `phaseinfo` data structure. 

**Note: Decompressing requires loading the original haplotype reference panel.** 


```julia
tgtfile = "mendel.imputed.jlso" # ultra-compressed genotypes after phasing & imputation
reffile = "ref.chr22.excludeTarget.vcf.gz" # original haplotype reference file
X1, X2, phaseinfo, sampleID, H = convert_compressed(Float64, tgtfile, reffile);
```

    [32mimporting reference data...100%|████████████████████████| Time: 0:01:52[39m


Check this compression protocol exhibit same error rate with [standard VCF compression](https://biona001.github.io/MendelImpute/dev/man/Phasing+and+Imputation/#Step-4:-%28only-for-simulated-data%29-check-imputation-accuracy). Note that `X1`, `X2`, and `H` are transposed. 


```julia
X_truth  = convert_gt(Float64, "target.chr22.full.vcf.gz") # import true genotypes
X_mendel = (X1 + X2)' # transpose X1 and X2
n, p = size(X_mendel)
println("error overall = $(sum(X_mendel .!= X_truth) / n / p)")
```

    error overall = 0.005397602533930585

