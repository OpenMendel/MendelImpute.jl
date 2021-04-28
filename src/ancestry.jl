###### This file is part of the MendelImpute.jl package.
###### It contains relevant code to perform chromosome painting and computing
###### computation of poplation admixtures.

"""
    admixture_global(tgtfile::String, reffile::String, 
        refID_to_population::Dict{String, String}, populations::Vector{String})

Computes global ancestry estimates for each sample in `tgtfile` using a labeled
reference panel `reffile`. 

# Inputs
- `tgtfile`: VCF or PLINK files. VCF files should end in `.vcf` or `.vcf.gz`.
    PLINK files should exclude `.bim/.bed/.fam` trailings but the trio must all
    be present in the same directory.
- `reffile`: Reference haplotype file ending in `.jlso` (compressed binary files).
    See [`compress_haplotypes`](@ref).
- `refID_to_population`: A dictionary mapping each sample IDs in the haplotype 
    reference panel to their population origin. For examples, see output of
    [`thousand_genome_population_to_superpopulation`](@ref) and
    [`thousand_genome_samples_to_super_population`](@ref)
- `populations`: A vector of `String` containing unique populations present in
    `values(refID_to_population)`. 

# Optional Inputs
- `Q_outfile`: Output file name for the estimated `Q` matrix. Default
    `Q_outfile="mendelimpute.ancestry.Q"`.
- `imputed_outfile`: Output file name for the imputed genotypes ending in `.jlso`.
    Default `impute_outfile = "mendelimpute.ancestry.Q.jlso"`

# Output
- `Q`: A `DataFrame` containing estimated ancestry fractions. Each row is a sample.
    Matrix will be saved in `mendelimpute.ancestry.Q`
"""
function admixture_global(
    tgtfile::AbstractString,
    reffile::AbstractString,
    refID_to_population::Dict{String, String},
    populations::Vector{String};
    Q_outfile = "mendelimpute.ancestry.Q",
    impute_outfile = "mendelimpute.ancestry.Q.jlso"
    )
    if !endswith(reffile, ".jlso")
        throw(ArgumentError("Reference file must be JLSO compressed!"))
    end

    # compute each person's phase information
    ph = phase(tgtfile, reffile, impute_outfile)

    # get ref sample IDs
    refID = MendelImpute.read_jlso(reffile).sampleID

    # compute sample composition using MendelImpute
    Q = DataFrame(zeros(length(ph), length(populations)), populations)
    for i in 1:length(ph)
        Q[i, :] .= composition(ph[i], refID, refID_to_population, populations=populations)
    end

    # save result in dataframe
    CSV.write(Q_outfile, Q)

    return Q
end

"""
    admixture_local(tgtfile::String, reffile::String, 
        refID_to_population::Dict{String, String}, populations::Vector{String},
        population_colors::Vector{RGB{FixedPointNumbers.N0f8}})

Computes global ancestry estimates for each sample in `tgtfile` using a labeled
reference panel `reffile`. 

# Inputs
- `tgtfile`: VCF or PLINK files. VCF files should end in `.vcf` or `.vcf.gz`.
    PLINK files should exclude `.bim/.bed/.fam` trailings but the trio must all
    be present in the same directory.
- `reffile`: Reference haplotype file ending in `.jlso` (compressed binary files).
    See [`compress_haplotypes`](@ref).
- `refID_to_population`: A dictionary mapping each sample IDs in the haplotype 
    reference panel to their population origin. For examples, see output of
    [`thousand_genome_population_to_superpopulation`](@ref) and
    [`thousand_genome_samples_to_super_population`](@ref)
- `population`: A list `String` containing unique populations present in
    `values(refID_to_population)`. 
- `population_colors`: A vector of colors for each population.
    `typeof(population_colors}` should be `Vector{RGB{FixedPointNumbers.N0f8}}`

# Output
- `Q`: Matrix containing estimated ancestry fractions. Each row is a haplotype.
    Sample 1's haplotypes are in rows 1 and 2, sample 2's are in rows 3, 4...etc.
- `pop_colors`: Matrix with sample dimension of `Q` storing colors. 
"""
function admixture_local(
    tgtfile::AbstractString,
    reffile::AbstractString,
    refID_to_population::Dict{String, String},
    populations::Vector{String},
    population_colors::Vector;
    outfile::AbstractString = "mendelimpute.ancestry.jlso"
    )
    if !endswith(reffile, ".jlso")
        throw(ArgumentError("Reference file must be JLSO compressed!"))
    end
    length(population_colors) == length(populations) || error("population_colors" * 
        "  should have same length as populations!")

    # compute each person's phase information
    ph = phase(tgtfile, reffile, outfile)

    # get ref sample IDs
    refID = MendelImpute.read_jlso(reffile).sampleID

    # compute each sample's composition
    sample_admixture = Vector{SampleAdmixture}(undef, length(ph))
    for i in 1:length(ph)
        sample_admixture[i] = paint(ph[i], refID, refID_to_population, populations=populations)
    end
    max_segments = maximum(length.(sample_admixture))

    # calculate admixture matrix and matrix of colors for easy plotting
    Q = zeros(2length(ph), max_segments)
    pop_colors = Matrix{eltype(population_colors)}(undef, 2length(ph), max_segments)
    for i in 1:length(ph)
        l1 = length(sample_admixture[i].H1)
        l2 = length(sample_admixture[i].H2)
        Q[2i - 1, 1:l1] .= sample_admixture[i].H1.composition
        Q[2i,     1:l2] .= sample_admixture[i].H2.composition
        for (j, pop) in enumerate(sample_admixture[i].H1.population_origins)
            pop_colors[2i - 1, j] = population_colors[findfirst(x -> x == pop, populations)]
        end
        for (j, pop) in enumerate(sample_admixture[i].H2.population_origins)
            pop_colors[2i, j] = population_colors[findfirst(x -> x == pop, populations)]
        end
    end

    return Q, pop_colors
end

"""
Data structure that stores one haplotype's local admixture information. The haplotype
is decomposed into segments whose individual length sums to 1. The lengths are stored 
in `composition`. Each segment `composition[i]` has an ancestry stored in 
`population_origins[i]`.
"""
struct HapolotypeAdmixture
    composition :: Vector{Float64}
    population_origins :: Vector{String}
end
Base.length(x::HapolotypeAdmixture) = length(x.composition)

"""
Data structure that stores a sample's (local) admixture information for both haplotypes
"""
struct SampleAdmixture
    H1 :: HapolotypeAdmixture
    H2 :: HapolotypeAdmixture
end

"""
Returns length of `x.H1` or `x.H2`, which ever is greater. 
"""
Base.length(x::SampleAdmixture) = max(length(x.H1), length(x.H2))

"""
    composition(sample_phase::HaplotypeMosaicPair, panelID::Vector{String}, 
        refID_to_population::Dict{String, String}, [populations::Vector{String}])

Computes a sample's chromosome composition based on phase information. This
function is used for easier plotting a person's admixed proportions.

# Inputs
- `sample_phase`: A `HaplotypeMosaicPair` storing phase information for a
    sample, includes haplotype start position and haplotype label.
- `panelID`: Sample ID's in the reference haplotype panel
- `refID_to_population`: A dictionary mapping each sample IDs in the haplotype 
    reference panel to their population origin. For examples, see output of
    [`thousand_genome_population_to_superpopulation`](@ref) and
    [`thousand_genome_samples_to_super_population`](@ref)

# Optional inputs
- `populations`: A unique list of populations present in `refID_to_population`

# Output
- `composition`: A list of percentages where `composition[i]` equals the
    sample's ancestry (in %) from `populations[i]` 
"""
function composition(
    sample_phase::HaplotypeMosaicPair,
    panelID::Vector{String},
    refID_to_population::Dict{String, String};
    populations::Vector{String} = unique(values(refID_to_population))
    )
    snps = sample_phase.strand1.length
    @assert snps == sample_phase.strand2.length "strands have different length!"

    # return elements
    composition = zeros(length(populations))

    # strand1 
    l1 = length(sample_phase.strand1.haplotypelabel)
    for i in 1:l1
        # get complete haplotype index and range of haplotype
        h1 = sample_phase.strand1.haplotypelabel[i]
        cur_range = (i == l1 ? 
            (sample_phase.strand1.start[i]:sample_phase.strand1.length) : 
            (sample_phase.strand1.start[i]:(sample_phase.strand1.start[i+1] 
            - 1)))

        # find current haplotype's population origin
        ref_id = panelID[div(h1, 2, RoundUp)]
        population = refID_to_population[ref_id]

        # update composition
        idx = findfirst(x -> x == population, populations)
        composition[idx] += length(cur_range)
    end

    # strand2
    l2 = length(sample_phase.strand2.haplotypelabel)
    for i in 1:l2
        # get complete haplotype index and range of haplotype
        h2 = sample_phase.strand2.haplotypelabel[i]
        cur_range = (i == l2 ? 
            (sample_phase.strand2.start[i]:sample_phase.strand2.length) : 
            (sample_phase.strand2.start[i]:(sample_phase.strand2.start[i+1] 
            - 1)))

        # find current haplotype's population origin
        ref_id = panelID[div(h2, 2, RoundUp)]
        population = refID_to_population[ref_id]

        # update composition
        idx = findfirst(x -> x == population, populations)
        composition[idx] += length(cur_range)
    end

    sum(composition) == 2snps || error("composition should sum to number of snps")
    
    return composition ./ 2snps
end

"""
    paint(sample_phase::HaplotypeMosaicPair, panelID::Vector{String},
        refID_to_population::Dict{String, String}, populations::Vector{String})

Converts a person's phased haplotype lengths into segments of percentages. This
function is used for easier plotting a "painted chromosome".

# Inputs
- `sample_phase`: A `HaplotypeMosaicPair` storing phase information for a
    sample, includes haplotype start position and haplotype label.
- `panelID`: Sample ID's in the reference haplotype panel
- `refID_to_population`: A dictionary mapping each sample IDs in the haplotype 
    reference panel to their population origin. For examples, see output of
    [`thousand_genome_population_to_superpopulation`](@ref) and
    [`thousand_genome_samples_to_super_population`](@ref)

# Optional inputs
- `populations`: A unique list of populations present in `refID_to_population`

# Output
- `s1_composition`: Haplotype 1's composition is stored in `s1_composition[1]`,
    which is a list of percentages such that `sum(s1_composition[1]) == 1`. Each
    segment `s1_composition[1][i]` has ancestry stored in `s1_composition[2][i]`.
- `s1_composition`: Haplotype 2's composition is stored in `s2_composition[1]`,
    which is a list of percentages such that `sum(s2_composition[1]) == 1`. Each
    segment `s2_composition[1][i]` has ancestry stored in `s2_composition[2][i]`.
"""
function paint(
    sample_phase::HaplotypeMosaicPair,
    panelID::Vector{String},
    refID_to_population::Dict{String, String};
    populations::Vector{String} = unique(values(refID_to_population))
    )
    snps = sample_phase.strand1.length
    @assert snps == sample_phase.strand2.length "strands have different length!"

    # length of strand 1 and 2
    l1 = length(sample_phase.strand1.haplotypelabel)
    l2 = length(sample_phase.strand2.haplotypelabel)

    # return elements
    s1_composition = (zeros(l1), Vector{String}(undef, l1))
    s2_composition = (zeros(l2), Vector{String}(undef, l2))

    # strand1 
    for i in 1:l1
        # get complete haplotype index and range of haplotype
        h1 = sample_phase.strand1.haplotypelabel[i]
        cur_range = (i == l1 ? 
            (sample_phase.strand1.start[i]:sample_phase.strand1.length) : 
            (sample_phase.strand1.start[i]:(sample_phase.strand1.start[i+1] 
            - 1)))

        # find current haplotype's population origin
        ref_id = panelID[div(h1, 2, RoundUp)]
        population = refID_to_population[ref_id]

        # update composition
        s1_composition[1][i] = length(cur_range) / snps
        s1_composition[2][i] = population
    end

    # strand2
    for i in 1:l2
        # get complete haplotype index and range of haplotype
        h2 = sample_phase.strand2.haplotypelabel[i]
        cur_range = (i == l2 ? 
            (sample_phase.strand2.start[i]:sample_phase.strand2.length) : 
            (sample_phase.strand2.start[i]:(sample_phase.strand2.start[i+1] 
            - 1)))

        # find current haplotype's population origin
        ref_id = panelID[div(h2, 2, RoundUp)]
        population = refID_to_population[ref_id]

        # update composition
        s2_composition[1][i] = length(cur_range) / snps
        s2_composition[2][i] = population
    end

    sum(s1_composition[1]) ≈ 1.0 || error("strand 1 compositions should sum to 1")
    sum(s2_composition[1]) ≈ 1.0 || error("strand 1 compositions should sum to 1")

    s1 = HapolotypeAdmixture(s1_composition[1], s1_composition[2])
    s2 = HapolotypeAdmixture(s2_composition[1], s2_composition[2])
    return SampleAdmixture(s1, s2)
    # return s1_composition, s2_composition
end

"""
    thousand_genome_samples_to_population()

Creates a dictionaries mapping sample IDs of 1000 genome project to 26 population
codes.

Population code and super population codes are described here:
https://www.internationalgenome.org/category/population/
"""
function thousand_genome_samples_to_population()
    # read population origin into a dataframe
    file = joinpath(normpath(MendelImpute.datadir()), "1000genomes.population.txt")
    df = CSV.read(file, DataFrame)

    # create dictionary with key = ID, value = population 
    refID_to_population = Dict{String, String}()
    for (id, population) in eachrow(df)
        refID_to_population[id] = population
    end
    refID_to_population
end

"""
    thousand_genome_samples_to_population()

Creates a dictionaries mapping sample IDs of 1000 genome project to 5
super population codes.

Population code and super population codes are described here:
https://www.internationalgenome.org/category/population/
"""
function thousand_genome_samples_to_super_population()
    # read population origin into a dataframe
    file = joinpath(normpath(MendelImpute.datadir()), "1000genomes.population.txt")
    df = CSV.read(file, DataFrame)

    # dict mapping population code to super population code
    pop_to_superpop = thousand_genome_population_to_superpopulation()

    # create dictionary with key = ID, value = population 
    refID_to_superpopulation = Dict{String, String}()
    for (id, population) in eachrow(df)
        refID_to_superpopulation[id] = pop_to_superpop[population]
    end
    refID_to_superpopulation
end

"""
    thousand_genome_population_to_superpopulation()

Creates a dictionary mapping population codes of 1000 genome project to their
super-population codes.

Population code and super population codes are described here:
https://www.internationalgenome.org/category/population/
"""
function thousand_genome_population_to_superpopulation()
    pop_to_superpop = Dict{String, String}()

    # 5 east asian
    pop_to_superpop["CHB"] = "EAS"
    pop_to_superpop["JPT"] = "EAS"
    pop_to_superpop["CHS"] = "EAS"
    pop_to_superpop["CDX"] = "EAS"
    pop_to_superpop["KHV"] = "EAS"

    # 5 european
    pop_to_superpop["CEU"] = "EUR"
    pop_to_superpop["TSI"] = "EUR"
    pop_to_superpop["FIN"] = "EUR"
    pop_to_superpop["GBR"] = "EUR"
    pop_to_superpop["IBS"] = "EUR"

    # 7 african
    pop_to_superpop["YRI"] = "AFR"
    pop_to_superpop["LWK"] = "AFR"
    pop_to_superpop["GWD"] = "AFR"
    pop_to_superpop["MSL"] = "AFR"
    pop_to_superpop["ESN"] = "AFR"
    pop_to_superpop["ASW"] = "AFR"
    pop_to_superpop["ACB"] = "AFR"

    # 4 ad mixed americans
    pop_to_superpop["MXL"] = "AMR"
    pop_to_superpop["PUR"] = "AMR"
    pop_to_superpop["CLM"] = "AMR"
    pop_to_superpop["PEL"] = "AMR"

    # 5 south asian
    pop_to_superpop["GIH"] = "SAS"
    pop_to_superpop["PJL"] = "SAS"
    pop_to_superpop["BEB"] = "SAS"
    pop_to_superpop["STU"] = "SAS"
    pop_to_superpop["ITU"] = "SAS"

    return pop_to_superpop
end
