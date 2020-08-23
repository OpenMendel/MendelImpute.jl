###### This file is part of the MendelImpute.jl package.
###### It contains relevant code to perform chromosome painting

# function paint(
#     sample_phase::HaplotypeMosaicPair,

#     )

# end

function paint(
    sample_phase::HaplotypeMosaicPair,
    sample_to_population::Dict{Symbol, Symbol},
    compressed_Hunique::CompressedHaplotypes;
    populations::Vector{Symbol} = unique_populations(sample_to_population),
    composition::AbstractVector = zeros(length(populations))
    )
    snps = sample_phase.strand1.length
    @assert snps == sample_phase.strand2.length "strands have different length!"
    
    l1 = length(sample_phase.strand1.haplotypelabel)
    l2 = length(sample_phase.strand2.haplotypelabel)
    ref_sampleIDs = compressed_Hunique.sampleID

    # strand1 
    for i in 1:l1
        h1 = sample_phase.strand1.haplotypelabel[i] # unique haplotype index
        w1 = sample_phase.strand1.window[i]
        cur_range = (i == l1 ? 
            (sample_phase.strand1.start[i]:sample_phase.strand1.length) : 
            (sample_phase.strand1.start[i]:(sample_phase.strand1.start[i+1] 
            - 1)))

        # convert unique haplotype idx to sampleID
        H1 = unique_all_idx_to_complete_idx(h1, w1, 
            compressed_Hunique) # complete haplotype idx
        sample_idx = div(H1, 2, RoundUp)
        id = Symbol(ref_sampleIDs[sample_idx])
        population = sample_to_population[id]

        # update composition
        idx = findfirst(x -> x == population, populations)
        composition[idx] += length(cur_range)
    end

    # strand2
    for i in 1:l2
        h2 = sample_phase.strand2.haplotypelabel[i] # unique haplotype index
        w2 = sample_phase.strand2.window[i]
        cur_range = (i == l2 ? 
            (sample_phase.strand2.start[i]:sample_phase.strand2.length) : 
            (sample_phase.strand2.start[i]:(sample_phase.strand2.start[i+1] 
            - 1)))

        # convert unique haplotype idx to sampleID
        H2 = unique_all_idx_to_complete_idx(h2, w2, 
            compressed_Hunique) # complete haplotype idx
        sample_idx = div(H2, 2, RoundUp)
        id = Symbol(ref_sampleIDs[sample_idx])
        population = sample_to_population[id]

        # update composition
        idx = findfirst(x -> x == population, populations)
        composition[idx] += length(cur_range)
    end

    # This is not strictly enforced since a sample's phase could have overlapping
    # regions due to breakpoint searching, which cases a very small region of snps
    # to be double counted.
    # sum(composition) == snps || error("composition should sum to number of snps")
    
    return composition ./ 2snps
end

"""
    unique_populations(x::Dict{Symbol, Symbol})

Computes the unique list of populations. `x` is a `Dict`
where each sample is a key and populations are values. 
"""
function unique_populations(x::Dict{Symbol, Symbol})
    populations = Symbol[]
    for (key, val) in x
        val ∉ populations && push!(populations, val)
    end
    return populations
end
