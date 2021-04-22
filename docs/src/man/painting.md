
# Estimating ancestry

If samples in the reference haplotype panel are labeled with a population origin, MendelImpute can also be used for:

+ Estimate admixed proportions
+ Chromosome painting


```julia
# first load all necessary packages
using Revise
using MendelImpute
using VCFTools
using GeneticVariation
using Random
using DataFrames
using Plots
using JLSO
using CSV
```

    ┌ Info: Precompiling MendelImpute [e47305d1-6a61-5370-bc5d-77554d143183]
    └ @ Base loading.jl:1317
    ┌ Info: Precompiling Plots [91a5bcdd-55d7-5caf-9e0b-520d859cae80]
    └ @ Base loading.jl:1317


## Prepare Example data for illustration

### Step 1. Filter chromosome data 

We use the [1000 genomes chromosome 22](http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/) as illustration.  The original data is filtered into target and reference panels. Follow [detailed example](https://openmendel.github.io/MendelImpute.jl/dev/man/Phasing+and+Imputation/#Detailed-Example) in Phasing and Imputation to obtain the same data.


!!! note

    In practice, it is better to infer ancestry of admixed populations using non-admixed reference populations. The example here is a simplified illustration and should not be taken too literally. 



### Step 2. Process each sample's population origin

MendelImpute needs to know each reference sample's origin (country/ethnicity/region...etc). This origin information should be provided by the reference haplotype panel, but users are free to further organize origin labels base on their own criteria. `MendelImpute` need a `Dict{key, value}` where each key is a reference sample ID and the value is the population code. Example dictionaries for 1000 genome project can be created by `MendelImpute`'s internal helper functions. Users not using 1000 genomes would have to manually construct such a dictionary mapping reference sample IDs to a desired population label. 

Here is a dictionary mapping reference sample IDs to super population codes. 


```julia
refID_to_superpopulation = thousand_genome_samples_to_super_population()
```




    Dict{String, String} with 2504 entries:
      "HG01791" => "EUR"
      "HG02736" => "SAS"
      "HG00182" => "EUR"
      "HG03914" => "SAS"
      "HG00149" => "EUR"
      "NA12156" => "EUR"
      "HG02642" => "AFR"
      "HG02851" => "AFR"
      "NA19835" => "AFR"
      "NA19019" => "AFR"
      "HG01131" => "AMR"
      "HG03578" => "AFR"
      "NA18550" => "EAS"
      "HG02401" => "EAS"
      "HG01350" => "AMR"
      "HG03973" => "SAS"
      "NA07000" => "EUR"
      "HG01709" => "EUR"
      "HG01395" => "AMR"
      "HG01980" => "AMR"
      "HG01979" => "AMR"
      "HG01122" => "AMR"
      "HG03869" => "SAS"
      "HG03729" => "SAS"
      "NA19920" => "AFR"
      ⋮         => ⋮



Here are dictionaries converting population code to super population codes.


```julia
pop_to_superpop = thousand_genome_population_to_superpopulation()
```




    Dict{String, String} with 26 entries:
      "CHS" => "EAS"
      "CDX" => "EAS"
      "GIH" => "SAS"
      "MSL" => "AFR"
      "KHV" => "EAS"
      "PUR" => "AMR"
      "ACB" => "AFR"
      "CLM" => "AMR"
      "FIN" => "EUR"
      "TSI" => "EUR"
      "BEB" => "SAS"
      "LWK" => "AFR"
      "STU" => "SAS"
      "JPT" => "EAS"
      "PJL" => "SAS"
      "ITU" => "SAS"
      "MXL" => "AMR"
      "GWD" => "AFR"
      "CEU" => "EUR"
      "YRI" => "AFR"
      "ASW" => "AFR"
      "ESN" => "AFR"
      "CHB" => "EAS"
      "IBS" => "EUR"
      "PEL" => "AMR"
      "GBR" => "EUR"



Note the [population codes](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/) for 1000 genome's samples are explained [here](https://www.internationalgenome.org/category/population/). 

## Estimate admixture proportions

+ The [composition](https://openmendel.github.io/MendelImpute.jl/dev/man/api/#MendelImpute.composition) will compute a list of percentages where `composition[i]` equals the sample's ancestry (in %) from `populations[i]`.
+ This illustration depends on **data preparation** above. 

### Step 1. Compute phase information using MendelImpute

This is equivalent to running a typical imputation. Please ensure that:
+ The output file name ends with `.jlso` (save output to ultra-compressed format)
+ `impute = true` (so the output contains the entire chromosome)

Note data used here is prepared in [Detailed Example](https://openmendel.github.io/MendelImpute.jl/dev/man/Phasing+and+Imputation/#Detailed-Example).


```julia
# compute each person's phase information
tgtfile = "target.chr22.typedOnly.masked.vcf.gz"
reffile = "ref.chr22.maxd1000.excludeTarget.jlso"
outfile = "mendel.imputed.jlso"
@time ph = phase(tgtfile, reffile, outfile);
```

    Number of threads = 1
    Importing reference haplotype data...


    [32mComputing optimal haplotypes...100%|████████████████████| Time: 0:00:26[39m
    [32mPhasing...100%|█████████████████████████████████████████| Time: 0:00:05[39m


    Total windows = 1634, averaging ~ 508 unique haplotypes per window.
    
    Timings: 
        Data import                     = 13.2588 seconds
            import target data             = 4.16205 seconds
            import compressed haplotypes   = 9.09677 seconds
        Computing haplotype pair        = 26.3008 seconds
            BLAS3 mul! to get M and N      = 1.10473 seconds per thread
            haplopair search               = 20.4958 seconds per thread
            initializing missing           = 0.1148 seconds per thread
            allocating and viewing         = 0.208966 seconds per thread
            index conversion               = 0.00789653 seconds per thread
        Phasing by win-win intersection = 5.64613 seconds
            Window-by-window intersection  = 0.571368 seconds per thread
            Breakpoint search              = 3.71102 seconds per thread
            Recording result               = 0.201057 seconds per thread
        Imputation                     = 3.85517 seconds
            Imputing missing               = 0.0245544 seconds
            Writing to file                = 3.83061 seconds
    
        Total time                      = 49.2067 seconds
    
     59.741822 seconds (118.80 M allocations: 7.179 GiB, 5.28% gc time)


### Step 2: import phased genotypes from JLSO file


```julia
# First import compressed reference panel
reffile = "ref.chr22.maxd1000.excludeTarget.jlso"
compressed_Hunique = MendelImpute.read_jlso(reffile)
panelID = compressed_Hunique.sampleID

# also need target sample's ancestry
tgtfile = "target.chr22.typedOnly.masked.vcf.gz"
reader = VCF.Reader(openvcf(tgtfile, "r"))
tgtID  = VCF.header(reader).sampleID
sample_superpopulation = [refID_to_superpopulation[id] for id in tgtID];
```


```julia
# here is each sample's super-population (sample 1 is EUR, sample 3 is EAS...etc)
sample_superpopulation
```




    100-element Vector{String}:
     "EUR"
     "EUR"
     "EAS"
     "EAS"
     "EAS"
     "EAS"
     "AMR"
     "AMR"
     "AMR"
     "AMR"
     "EUR"
     "AMR"
     "EUR"
     ⋮
     "AMR"
     "AFR"
     "AFR"
     "EUR"
     "EUR"
     "EUR"
     "EUR"
     "EUR"
     "EUR"
     "EUR"
     "SAS"
     "SAS"



### Step 3: call `composition` function

+ The [composition](https://openmendel.github.io/MendelImpute.jl/dev/man/api/#MendelImpute.composition) will compute a list of percentages where `composition[i]` equals the sample's ancestry (in %) from `populations[i]`.
+ We are finally using the imputation result stored in `ph`.


```julia
populations = MendelImpute.unique_populations(refID_to_superpopulation)
@time sample1_comp = composition(ph[1], panelID, refID_to_superpopulation) # origin GBR (EUR)
@time sample4_comp = composition(ph[4], panelID, refID_to_superpopulation) # origin CHS (EAS)
@time sample84_comp = composition(ph[84], panelID, refID_to_superpopulation) # origin LWK (AFR)

println("sample 1 = ", round(sample1_comp[1], digits=3), " S. asian")
println("sample 1 = ", round(sample1_comp[2], digits=3), " E. asian")
println("sample 1 = ", round(sample1_comp[3], digits=3), " European")
println("sample 1 = ", round(sample1_comp[4], digits=3), " Admixed-American")
println("sample 1 = ", round(sample1_comp[5], digits=3), " Africans\n")

println("sample 4 = ", round(sample4_comp[1], digits=3), " S. asian")
println("sample 4 = ", round(sample4_comp[2], digits=3), " E. asian")
println("sample 4 = ", round(sample4_comp[3], digits=3), " European")
println("sample 4 = ", round(sample4_comp[4], digits=3), " Admixed-American")
println("sample 4 = ", round(sample4_comp[5], digits=3), " Africans\n")
    
println("sample 84 = ", round(sample84_comp[1], digits=3), " S. asian")
println("sample 84 = ", round(sample84_comp[2], digits=3), " E. asian")
println("sample 84 = ", round(sample84_comp[3], digits=3), " European")
println("sample 84 = ", round(sample84_comp[4], digits=3), " Admixed-American")
println("sample 84 = ", round(sample84_comp[5], digits=3), " Africans");
```

      0.004967 seconds (24 allocations: 1.703 KiB, 96.70% compilation time)
      0.000163 seconds (6 allocations: 544 bytes)
      0.000177 seconds (6 allocations: 544 bytes)
    sample 1 = 0.652 S. asian
    sample 1 = 0.089 E. asian
    sample 1 = 0.023 European
    sample 1 = 0.169 Admixed-American
    sample 1 = 0.068 Africans
    
    sample 4 = 0.189 S. asian
    sample 4 = 0.061 E. asian
    sample 4 = 0.01 European
    sample 4 = 0.053 Admixed-American
    sample 4 = 0.687 Africans
    
    sample 84 = 0.065 S. asian
    sample 84 = 0.014 E. asian
    sample 84 = 0.784 European
    sample 84 = 0.111 Admixed-American
    sample 84 = 0.025 Africans


Here `sample1_comp[i]` equals the sample's estimated ancestry (in %) from `populations[i]`. 

**Conclusion:** We computed the population percentages for sample 1, 4, and 84 with respect to the 5 reference super populations. Thus sample 1 is 65% European, 10% South Asian, 20% American...etc. Sample 4 is 20% European, 70% East Asian,...etc. Sample 84 is 80% African and 5% European...etc. 

## Chromosome painting

The main function is the [paint](https://openmendel.github.io/MendelImpute.jl/dev/man/api/#MendelImpute.paint) function. For an imputed sample, it will convert **each haplotype segment** into a percentage indicating the segment's length in the chromosome. Then the list can be used for easy plotting. 

**Note:** this illustration depends on **data preparation** above. 

### Step 1: Choose your colors

In this example, colors are arranged such that:
+ Blue ≈ European/American
+ Red ≈ South/East Asian
+ Green ≈ African

Of course, Julia lets you plot your favoriate colors. We pick our colors here: https://mdigi.tools/color-shades/#008000.


```julia
continent = ["SAS", "EAS", "EUR", "AMR", "AFR"]
continent_colors = [colorant"#e6194B", colorant"#800000", colorant"#4363d8", colorant"#0000b3", colorant"#bfef45"]
```




![svg](output_17_0.svg)



### Step 2: Run `paint` funcion

This function convert the imputed haplotype segments into a list of percentages (one list for each strand). This is simply a post-processing routine so that data can be used for easy plotting later.


```julia
populations = unique_populations(refID_to_superpopulation)
@time sample1_s1_comp, sample1_s2_comp = paint(ph[1], panelID, refID_to_superpopulation, populations=populations)
@time sample4_s1_comp, sample4_s2_comp = paint(ph[4], panelID, refID_to_superpopulation, populations=populations)
@time sample84_s1_comp, sample84_s2_comp = paint(ph[84], panelID, refID_to_superpopulation, populations=populations);
```

      0.080868 seconds (114.69 k allocations: 6.641 MiB, 98.85% compilation time)
      0.000094 seconds (12 allocations: 19.906 KiB)
      0.000101 seconds (12 allocations: 22.406 KiB)



```julia
# view phase information of haplotype 1 from sample 1
[sample1_s1_comp[1] sample1_s1_comp[2]]
```




    545×2 Matrix{Any}:
     0.000562617  "AFR"
     0.000447699  "AMR"
     0.000476429  "EUR"
     0.0002849    "EUR"
     0.0001221    "EUR"
     7.66117e-5   "EUR"
     0.000287294  "EUR"
     0.000411788  "EAS"
     0.000205894  "EAS"
     0.00117072   "EUR"
     0.00320811   "EUR"
     0.00117312   "EUR"
     0.00103426   "EUR"
     ⋮            
     0.00115636   "AMR"
     0.00370609   "EUR"
     0.000323205  "EAS"
     0.00120424   "EUR"
     0.000433335  "EUR"
     0.00104623   "EUR"
     0.00195839   "EUR"
     0.00037827   "EUR"
     0.000926522  "EUR"
     0.000938493  "EAS"
     0.000646411  "EUR"
     0.00162321   "EUR"



### Step 3: Generate plots for painted chromosomes

We found the [StatsPlots.jl](https://github.com/JuliaPlots/StatsPlots.jl) (version 0.14.17) package to be more useful for this purpose, although the code below still did the plotting in a very roundabout way. Newer StatsPlots version breaks this plotting code. To pin `StatsPlots` to the correct version, execute `]pin StatsPlots@v0.14.17` in Julia.


```julia
using StatsPlots, FixedPointNumbers

# assign a color to each haplotype segment
sample1_s1_colors = [continent_colors[findfirst(x -> x == pop, continent)] for pop in sample1_s1_comp[2]]
sample1_s1_colors = reshape(sample1_s1_colors, 1, length(sample1_s1_colors))
sample1_s2_colors = [continent_colors[findfirst(x -> x == pop, continent)] for pop in sample1_s2_comp[2]]
sample1_s2_colors = reshape(sample1_s2_colors, 1, length(sample1_s2_colors))
sample4_s1_colors = [continent_colors[findfirst(x -> x == pop, continent)] for pop in sample4_s1_comp[2]]
sample4_s1_colors = reshape(sample4_s1_colors, 1, length(sample4_s1_colors))
sample4_s2_colors = [continent_colors[findfirst(x -> x == pop, continent)] for pop in sample4_s2_comp[2]]
sample4_s2_colors = reshape(sample4_s2_colors, 1, length(sample4_s2_colors))
sample84_s1_colors = [continent_colors[findfirst(x -> x == pop, continent)] for pop in sample84_s1_comp[2]]
sample84_s1_colors = reshape(sample84_s1_colors, 1, length(sample84_s1_colors))
sample84_s2_colors = [continent_colors[findfirst(x -> x == pop, continent)] for pop in sample84_s2_comp[2]]
sample84_s2_colors = reshape(sample84_s2_colors, 1, length(sample84_s2_colors));

# roundabout code for plotting...
sample1_s1l = length(sample1_s1_comp[1])
sample1_s2l = length(sample1_s2_comp[1])
sample4_s1l = length(sample4_s1_comp[1])
sample4_s2l = length(sample4_s2_comp[1])
sample84_s1l = length(sample84_s1_comp[1])
sample84_s2l = length(sample84_s2_comp[1])
maxlen = max(sample1_s1l, sample1_s2l, sample4_s1l, sample4_s2l, sample84_s1l, sample84_s2l)

mydata = zeros(6, maxlen)
copyto!(@view(mydata[1, 1:sample1_s1l]), sample1_s1_comp[1])
copyto!(@view(mydata[2, 1:sample1_s2l]), sample1_s2_comp[1])
copyto!(@view(mydata[3, 1:sample4_s1l]), sample4_s1_comp[1])
copyto!(@view(mydata[4, 1:sample4_s2l]), sample4_s2_comp[1])
copyto!(@view(mydata[5, 1:sample84_s1l]), sample84_s1_comp[1])
copyto!(@view(mydata[6, 1:sample84_s2l]), sample84_s2_comp[1])

pop_colors = Matrix{RGB{Normed{UInt8,8}}}(undef, 6, maxlen)
copyto!(@view(pop_colors[1, 1:sample1_s1l]), sample1_s1_colors)
copyto!(@view(pop_colors[2, 1:sample1_s2l]), sample1_s2_colors)
copyto!(@view(pop_colors[3, 1:sample4_s1l]), sample4_s1_colors)
copyto!(@view(pop_colors[4, 1:sample4_s2l]), sample4_s2_colors)
copyto!(@view(pop_colors[5, 1:sample84_s1l]), sample84_s1_colors)
copyto!(@view(pop_colors[6, 1:sample84_s2l]), sample84_s2_colors)

xnames = ["Sample 1 hap1", "Sample 1 hap2", "Sample 4 hap1", "Sample 4 hap2", "Sample 84 hap1", "Sample 84 hap2"]
ynames = ["SNP 1", "SNP 208k", "SNP 417k"]

# color haplotypes
chrom_plt2 = groupedbar(mydata, bar_position = :stack, bar_width=0.7, label=:none, 
    lw = 0, color=pop_colors, xticks=(1:1:6, xnames), yticks=(0:0.5:1, ynames),
    ytickfont=font(12), xtickfont=font(12), xrotation=20, right_margin = 35Plots.mm,
    grid=false)

# create a plot for legend
color_x = ones(5)
color_y = collect(1:1:5)
admixture_chrom_plt = scatter!(color_x, color_y, color=reverse(continent_colors), ytick=(1:1:5, reverse(continent)), 
    xrange=(0.9, 1.1), xtick=false, label=:none, markersize=8, ytickfont=font(16),
    grid=false, framestyle=:grid, mirror=true, tick_direction=:out, markershape=:rect,
    inset = (1, bbox(-0.05, 0.0, 0.05, 1.0, :bottom, :right)), subplot = 2)
```




![svg](output_22_0.svg)




```julia
]st StatsPlots # ensure StatsPlots version is v0.14.17
```

    [36m[1m     Project[22m[39m MendelImpute v1.1.0
    [32m[1m      Status[22m[39m `~/.julia/dev/MendelImpute/Project.toml`
     [90m [f3b207a7] [39m[37mStatsPlots v0.14.17 ⚲[39m


**Conclusion:** 
+ We can visualize the linkage patterns for the 3 samples across their 6 haplotypes
+ Sample 1 is mostly European and admixed American, sample 2 is mainly South/East Asian, and sample 3 is mainly African.

**Note:** this example should not be taken too literally, since we *did not* exclude admixed samples from the reference panel. For more details, please refer to our paper, or file an issue on GitHub. 
