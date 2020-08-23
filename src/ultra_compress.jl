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
    loaded = JLSO.load(phaseinfo)
    phase = loaded[:ph]

    X1, X2 = convert(t, phase, H)
    return X1, X2, phase[:sampleID]
end

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

"""
Same as update_phase! but records complete haplotype index. Won't record
if previous haplotype label is the same as current one. This helper function
may not be needed. 
"""
function update_compressed_phase!(ph::HaplotypeMosaic, bkpt::Int, hap_prev,
    hap_curr, Xwi_start::Int, Xwi_mid::Int, Xwi_end::Int)

    X_bkpt_end = Xwi_start + bkpt

    # no breakpoints or double breakpoints
    if bkpt == -1
        ph.haplotypelabel[end] == hap_curr && return nothing # only push new segment
        push!(ph.start, Xwi_mid)
        push!(ph.haplotypelabel, hap_curr)
        return nothing
    end

    # previous window's haplotype completely covers current window
    if bkpt == length(Xwi_start:Xwi_end)
        ph.haplotypelabel[end] == hap_prev && return nothing # only push new segment
        push!(ph.start, Xwi_mid)
        push!(ph.haplotypelabel, hap_prev)
        return nothing
    end

    if Xwi_mid <= X_bkpt_end <= Xwi_end
        # previous window extends to current window
        if ph.haplotypelabel[end] != hap_prev
            push!(ph.start, Xwi_mid)
            push!(ph.haplotypelabel, hap_prev)
        end
        # 2nd part of current window
        if ph.haplotypelabel[end] != hap_curr
            push!(ph.start, X_bkpt_end)
            push!(ph.haplotypelabel, hap_curr)
        end
    elseif X_bkpt_end < Xwi_mid
        # current window extends to previous window
        if ph.haplotypelabel[end] != hap_curr
            push!(ph.start, X_bkpt_end)
            push!(ph.haplotypelabel, hap_curr)
        end
        # update current window
        if ph.haplotypelabel[end] != hap_curr
            push!(ph.start, Xwi_mid)
            push!(ph.haplotypelabel, hap_curr)
        end
    else
        error("bkpt does not satisfy -1 <= bkpt <= 2width!")
    end

    return nothing
end