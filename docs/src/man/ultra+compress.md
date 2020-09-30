
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
@time phaseinfo = phase(tgtfile, reffile, outfile);
```

    Number of threads = 1
    Importing reference haplotype data...


    [32mComputing optimal haplotypes...100%|████████████████████| Time: 0:00:23[39m
    [32mPhasing...100%|█████████████████████████████████████████| Time: 0:00:05[39m


    Total windows = 1634, averaging ~ 508 unique haplotypes per window.
    
    Timings: 
        Data import                     = 13.5855 seconds
            import target data             = 2.99557 seconds
            import compressed haplotypes   = 10.5899 seconds
        Computing haplotype pair        = 23.5138 seconds
            BLAS3 mul! to get M and N      = 1.02228 seconds per thread
            haplopair search               = 18.3673 seconds per thread
            initializing missing           = 0.101495 seconds per thread
            allocating and viewing         = 0.270367 seconds per thread
            index conversion               = 0.0212053 seconds per thread
        Phasing by win-win intersection = 5.17088 seconds
            Window-by-window intersection  = 0.499549 seconds per thread
            Breakpoint search              = 3.71212 seconds per thread
            Recording result               = 0.00855205 seconds per thread
        Imputation                     = 3.04347 seconds
            Imputing missing               = 0.0513904 seconds
            Writing to file                = 2.99208 seconds
    
        Total time                      = 45.5041 seconds
    
     60.231372 seconds (124.56 M allocations: 6.857 GiB, 5.31% gc time)


The object saved to `mendel.imputed.jlso` is literally the `phaseinfo` variable. We can inspect its element:


```julia
# look at sample 1's haplotype segments
haplotype_labels = phaseinfo[1].strand1.haplotypelabel # strand1
haplotype_start = phaseinfo[1].strand1.start # strand1
[haplotype_start haplotype_labels]
```




    545×2 Array{Int64,2}:
          1  4119
        236   887
        423   272
        622    12
        741   124
        792     4
        824    24
        944  1282
       1116  1741
       1202  4543
       1691  1198
       3031    22
       3521    18
          ⋮  
     411702   877
     412185    74
     413733  3849
     413868   248
     414371    31
     414552  3187
     414989  4481
     415807     5
     415965   143
     416352  1276
     416744    71
     417014   311



## Loading

The function [convert_compressed](https://OpenMendel.github.io/MendelImpute.jl/dev/man/api/#MendelImpute.convert_compressed) will load the ultra-compressed data into genotype matrices and the original `phaseinfo` data structure. 

**Note: Decompressing requires loading the original haplotype reference panel.** 


```julia
tgtfile = "mendel.imputed.jlso" # ultra-compressed genotypes after phasing & imputation
reffile = "ref.chr22.excludeTarget.vcf.gz" # original haplotype reference file
X1, X2, phaseinfo, sampleID, H = convert_compressed(Float64, tgtfile, reffile);
```

    [32mimporting reference data...100%|████████████████████████| Time: 0:01:51[39m


Check this compression protocol exhibit same error rate with [standard VCF compression](https://OpenMendel.github.io/MendelImpute.jl/dev/man/Phasing+and+Imputation/#Step-4:-%28only-for-simulated-data%29-check-imputation-accuracy). Note that `X1`, `X2`, and `H` are transposed. 


```julia
X_truth  = convert_gt(Float64, "target.chr22.full.vcf.gz") # import true genotypes
X_mendel = (X1 + X2)' # transpose X1 and X2
n, p = size(X_mendel)
println("error overall = $(sum(X_mendel .!= X_truth) / n / p)")
```

    error overall = 0.00527504782243333

