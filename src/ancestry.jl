###### This file is part of the MendelImpute.jl package.
###### It contains relevant code to perform chromosome painting and computing
###### computation of poplation admixtures.

"""
    composition(sample_phase::HaplotypeMosaicPair, panelID::Vector{String}, 
        refID_to_population::Dict{String, String}, [populations::Vector{String}])

Computes a sample's chromosome composition based on phase information. This
function is used for easier plotting a person's admixed proportions.

# Inputs
- `sample_phase`: A `HaplotypeMosaicPair` storing phase information for a
    sample, includes haplotype start position and haplotype label.
- `panelID`: Sample ID's in the reference haplotype panel
- `refID_to_population`: A dictionary mapping each ID in the haplotype 
    reference panel to its population origin. 

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
    populations::Vector{String} = unique_populations(refID_to_population)
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
    unique_populations(x::Dict{String, String})

Computes the unique list of populations, preserving order. `x` is a `Dict`
where each sample is a key and populations are values. 
"""
function unique_populations(x::Dict{String, String})
    populations = String[]
    for (key, val) in x
        val ∉ populations && push!(populations, val)
    end
    return populations
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
- `refID_to_population`: A dictionary mapping each ID in the haplotype 
    reference panel to its population origin. 

# Optional inputs
- `populations`: A unique list of populations present in `refID_to_population`

# Output
- `composition`: A list of percentages where `composition[i]` equals the
    sample's ancestry (in %) from `populations[i]` 
"""
function paint(
    sample_phase::HaplotypeMosaicPair,
    panelID::Vector{String},
    refID_to_population::Dict{String, String};
    populations::Vector{String} = unique_populations(refID_to_population)
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

    return s1_composition, s2_composition
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
