###### This file is part of the MendelImpute.jl package.
###### It contains relevant code to decompress a vector of `HaplotypeMosaicPair`
###### to a genotype matrix

"""
convert_compressed(t, phaseinfo, H)

Converts `phaseinfo` into a phased genotype matrix of type `t` using the full
reference haplotype panel `H` 

# Inputs
- `t`: Type of matrix. If `bool`, genotypes are converted to a `BitMatrix`
- `phaseinfo`: Vector of `HaplotypeMosaicPair`s stored in `.jlso` format
- `reffile`: The complete (uncompressed) haplotype reference file

# Output
- `(X1, X2)`: Tuple of matrix where `X1` is allele1 and `X2` is allele2. Each 
    column is a sample. 
- `sampleID`: The ID's of each imputed person. 
"""
function convert_compressed(
    t::Type{T}, 
    phaseinfo::AbstractString,
    reffile::AbstractString
    ) where T <: Real
    endswith(phaseinfo, ".jlso") || error("phaseinfo does not end with '.jlso'")
    H = convert_ht(Bool, reffile, trans=true, msg="importing reference data...")
    
    phase_and_sampleID = JLSO.load(phaseinfo)
    phase = phase_and_sampleID[:ph]
    sampleID = phase_and_sampleID[:sampleID]

    X1, X2 = convert_compressed(t, phase, H)
    return X1, X2, phase, sampleID, H
end

"""
    convert_compressed(t, phaseinfo, H)

Columns of `H` are haplotypes.
"""
function convert_compressed(
    t::Type{T},
    phaseinfo::Vector{HaplotypeMosaicPair},
    H::AbstractMatrix
    ) where T <: Real

    people = length(phaseinfo)
    snps = phaseinfo[1].strand1.length
    M = (t == Bool ? BitArray{2} : Matrix{t})

    X1 = M(undef, snps, people)
    X2 = M(undef, snps, people)
    impute!(X1, X2, H, phaseinfo)

    return X1, X2
end

# function complete_idx_to_unique_all_idx!(
#     sample_phase::HaplotypeMosaicPair,
#     compressed_Hunique::CompressedHaplotypes,
#     )

#     l1 = length(sample_phase.strand1.haplotypelabel)
#     l2 = length(sample_phase.strand2.haplotypelabel)

#     # strand1 
#     for i in 1:l1
#         h1 = sample_phase.strand1.haplotypelabel[i] # unique haplotype index
#         w1 = sample_phase.strand1.window[i]
#         H1 = unique_all_idx_to_complete_idx(h1, w1, 
#             compressed_Hunique) # complete haplotype idx
#         sample_phase.strand1.haplotypelabel[i] = H1
#     end

#     # strand2
#     for i in 1:l2
#         h2 = sample_phase.strand2.haplotypelabel[i] # unique haplotype index
#         w2 = sample_phase.strand2.window[i]
#         H2 = unique_all_idx_to_complete_idx(h2, w2, 
#             compressed_Hunique) # complete haplotype idx 
#         sample_phase.strand2.haplotypelabel[i] = H2
#     end
# end
