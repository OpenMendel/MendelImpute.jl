function convert_gt(
    b::Bgen,
    T=Float64
    )
    n = n_samples(b)
    p = n_variants(b)
    
    # return arrays
    G = Matrix{T}(undef, p, n)
    sampleID = Vector{String}(undef, n)
    chr = Vector{String}(undef, p)
    pos = Vector{String}(undef, p)
    snpID = Vector{String}(undef, p)
    ref = Vector{String}(undef, p)
    alt = Vector{String}(undef, p)

    # loop over each variant
    i = 1
    for v in iterator(b; from_bgen_start=true)
        dose = minor_allele_dosage!(b, v; T=T)
        copyto!(@view(G[i, :]), dose)
        chr[i], pos[i], snpID[i], ref[i], alt[i] = chrom(v), pos(v), rsid(v),
            major_allele(v), minor_allele(v)
        i += 1
        clear!(v)
    end
    return G, sampleID, chr, pos, snpID, ref, alt
end
