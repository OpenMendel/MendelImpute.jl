###### This file is part of the MendelImpute.jl package.
###### It contains relevant code to perform chromosome painting as well as
###### code to decompress a vector of `HaplotypeMosaicPair` to a genotype matrix

function paint(
    sample_phase::HaplotypeMosaicPair,
    sample_to_population::Dict{Symbol, Symbol},
    compressed_Hunique::CompressedHaplotypes;
    populations::Vector{Symbol} = unique_populations(sample_to_population),
    composition::AbstractVector = zeros(length(populations))
    )
    snps = sample_phase.strand1.length
    l = length(sample_phase.strand1.haplotypelabel)
    ref_sampleIDs = compressed_Hunique.sampleID

    # strand1 
    for i in 1:l
        h1 = sample_phase.strand1.haplotypelabel[i] # unique haplotype index
        w1 = sample_phase.strand1.window[i]
        snps = (i == l ? 
            (sample_phase.strand1.start[i]:sample_phase.strand1.length) : 
            (sample_phase.strand1.start[i]:(sample_phase.strand1.start[i+1] 
            - 1)))

        # convert unique haplotype idx to sampleID
        H1 = unique_all_idx_to_complete_idx(h1, w1, 
            compressed_Hunique) # complete haplotype idx
        sample_idx = div(H1, 2, RoundUp)
        id = ref_sampleIDs[sample_idx]
        population = Symbol(sample_to_population[id])

        # update composition
        idx = findfirst(x -> x == population, populations)
        composition[idx] += length(snps)
    end

    # strand2
    for i in 1:l
        h2 = sample_phase.strand2.haplotypelabel[i] # unique haplotype index
        w2 = sample_phase.strand2.window[i]
        snps = (i == l ? 
            (sample_phase.strand2.start[i]:sample_phase.strand2.length) : 
            (sample_phase.strand2.start[i]:(sample_phase.strand2.start[i+1] 
            - 1)))

        # convert unique haplotype idx to sampleID
        H1 = unique_all_idx_to_complete_idx(h2, w2, 
            compressed_Hunique) # complete haplotype idx
        sample_idx = div(H1, 2, RoundUp)
        id = ref_sampleIDs[sample_idx]
        population = Symbol(sample_to_population[id])

        # update composition
        idx = findfirst(x -> x == population, populations)
        composition[idx] += length(snps)
    end
    
    return composition
end

"""
unique_populations(x::Dict{Symbol, Symbol})

Computes the unique list of populations. `x` is a `Dict`
where each sample is a key and populations are values. 
"""
function unique_populations(x::Dict{Symbol, Symbol})
    populations = Symbol[]
    for (key, val) in x
        if val ∉ populations 
            num_unique += 1
            push!(populations, val)
        end
    end
    return populations
end

# """
#     convert(t, phaseinfo, compressed_Hunique)

# Converts haplotype segments stored in `phaseinfo` with respect to 
# `compressed_Hunique` into a phased genotype matrix of type `t`

# # Inputs
# - `t`: Type of matrix. If `bool`, genotypes are converted to a `BitMatrix`
# - `phaseinfo`: Vector of `HaplotypeMosaicPair`s that store haplotype segments
# - `CompressedHaplotypes` A `.jlso` compressed haplotype reference panel for 
#     which `phaseinfo` is recorded with respect to. 

# # Output
# - `(X1, X2)`: Tuple of matrix where `X1` is allele1 and `X2` is allele2. 
# """
# function convert(
#     t::Type{T},
#     phaseinfo::Vector{HaplotypeMosaicPair},
#     compressed_Hunique::CompressedHaplotypes
#     ) where T <: Real
#     people = length(phaseinfo)
#     snps = phaseinfo[1].strand1.length
#     M = (t == Bool ? BitArray{2} : Matrix{t})

#     X1 = M(undef, snps, people)
#     X2 = M(undef, snps, people)
#     impute!(X1, X2, compressed_Hunique, phaseinfo)

#     return X1, X2
# end

# function convert(
#     t::Type{T}, 
#     phaseinfo::AbstractString,
#     compressed_Hunique::CompressedHaplotypes
#     ) where T <: Real
#     endswith(phaseinfo, ".jlso") || error("phaseinfo does not end with '.jlso'")
#     return convert(t, JLSO.load(phaseinfo)[:ph], compressed_Hunique)
# end

# function convert(
#     t::Type{T}, 
#     phaseinfo::AbstractString,
#     compressed_haplotypes::AbstractString
#     ) where T <: Real
#     endswith(phaseinfo, ".jlso") || error("phaseinfo does not end with '.jlso'")
#     endswith(compressed_haplotypes, ".jlso") || 
#         error("compressed_haplotypes does not end with '.jlso'")
#     return convert(t, JLSO.load(phaseinfo)[:ph], 
#         JLSO.load(compressed_haplotypes)[:compressed_Hunique])
# end
