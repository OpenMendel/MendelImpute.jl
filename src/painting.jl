###### This file is part of the MendelImpute.jl package.
###### It contains relevant code to perform chromosome painting as well as
###### code to decompress a vector of `HaplotypeMosaicPair` to a genotype matrix

"""
    convert(t, phaseinfo, compressed_Hunique)

Converts haplotype segments stored in `phaseinfo` with respect to 
`compressed_Hunique` into a phased genotype matrix of type `t`

# Inputs
- `t`: Type of matrix. If `bool`, genotypes are converted to a `BitMatrix`
- `phaseinfo`: Vector of `HaplotypeMosaicPair`s that store haplotype segments
- `CompressedHaplotypes` A `.jlso` compressed haplotype reference panel for 
    which `phaseinfo` is recorded with respect to. 

# Output
- `(X1, X2)`: Tuple of matrix where `X1` is allele1 and `X2` is allele2. 
"""
function convert(
    t::Type{T},
    phaseinfo::Vector{HaplotypeMosaicPair},
    compressed_Hunique::CompressedHaplotypes
    ) where T <: Real
    people = length(phaseinfo)
    snps = phaseinfo[1].strand1.length
    M = (t == Bool ? BitArray{2} : Matrix{t})

    X1 = M(undef, snps, people)
    X2 = M(undef, snps, people)
    impute!(X1, X2, compressed_Hunique, phaseinfo)

    return X1, X2
end

function convert(
    t::Type{T}, 
    phaseinfo::AbstractString,
    compressed_Hunique::CompressedHaplotypes
    ) where T <: Real
    endswith(phaseinfo, ".jlso") || error("phaseinfo does not end with '.jlso'")
    return convert(t, JLSO.load(phaseinfo)[:ph], compressed_Hunique)
end

function convert(
    t::Type{T}, 
    phaseinfo::AbstractString,
    compressed_haplotypes::AbstractString
    ) where T <: Real
    endswith(phaseinfo, ".jlso") || error("phaseinfo does not end with '.jlso'")
    endswith(compressed_haplotypes, ".jlso") || 
        error("compressed_haplotypes does not end with '.jlso'")
    return convert(t, JLSO.load(phaseinfo)[:ph], 
        JLSO.load(compressed_haplotypes)[:compressed_Hunique])
end
