"""
    impute!(X, compressed_haplotypes, phaseinfo, outfile, X_sampleID)

Imputes `X` using `phaseinfo` and outputs result in `outfile`. All genotypes 
in `outfile` are non-missing. If `XtoH_idx == nothing`, all SNPs in reference
file will be imputed. 
"""
function impute!(
    X::AbstractMatrix,
    compressed_haplotypes::CompressedHaplotypes,
    phaseinfo::Vector{HaplotypeMosaicPair},
    outfile::AbstractString,
    X_sampleID::AbstractVector;
    XtoH_idx::Union{Nothing, AbstractVector} = nothing
    )
    # impute without changing observed entries
    impute_discard_phase!(X, compressed_haplotypes, phaseinfo)

    # retrieve reference file information
    chr = (isnothing(XtoH_idx) ? compressed_haplotypes.chr : compressed_haplotypes.chr[XtoH_idx])
    pos = (isnothing(XtoH_idx) ? compressed_haplotypes.pos : compressed_haplotypes.pos[XtoH_idx])
    ids = (isnothing(XtoH_idx) ? compressed_haplotypes.SNPid : compressed_haplotypes.SNPid[XtoH_idx])
    ref = (isnothing(XtoH_idx) ? compressed_haplotypes.refallele : compressed_haplotypes.refallele[XtoH_idx])
    alt = (isnothing(XtoH_idx) ? compressed_haplotypes.altallele : compressed_haplotypes.altallele[XtoH_idx])

    # write minimal meta information to outfile
    io = openvcf(outfile, "w")
    pb = PipeBuffer()
    print(pb, "##fileformat=VCFv4.2\n")
    print(pb, "##source=MendelImpute\n")
    print(pb, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

    # header line should match reffile (i.e. sample ID's should match)
    print(pb, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    for id in X_sampleID
        print(pb, "\t", id)
    end
    print(pb, "\n")
    (bytesavailable(pb) > (16*1024)) && write(io, take!(pb))

    pmeter = Progress(size(X, 1), 5, "Writing to file...")
    @inbounds for i in 1:size(X, 1)
        # write meta info (chrom/pos/id/ref/alt)
        print(pb, chr[i], "\t", string(pos[i]), "\t", ids[i][1], "\t", ref[i], "\t", alt[i][1], "\t.\tPASS\t.\tGT")
        
        for j in 1:size(X, 2)
            if X[i, j] == 0
                print(pb, "\t0/0")
            elseif X[i, j] == 1
                print(pb, "\t1/0")
            elseif X[i, j] == 2
                print(pb, "\t1/1")
            else
                error("imputed genotypes can only be 0, 1, 2 but X[$i, $j]) = $(X[i, j])")
            end
        end
        print(pb, "\n")
        (bytesavailable(pb) > (16*1024)) && write(io, take!(pb))
        next!(pmeter)
    end
    write(io, take!(pb))

    # close & return
    close(io); close(pb)
    return nothing
end

"""
    impute!(X, H, phase)

Imputes `X` completely using segments of haplotypes `H` where segments are stored in `phase`. 
Non-missing entries in `X` can be different after imputation. Preserves phase information.
"""
function impute!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    phase::Vector{HaplotypeMosaicPair}
    )

    fill!(X, 0)
    # loop over individuals
    for i in 1:size(X, 2)
        for s in 1:(length(phase[i].strand1.start) - 1)
            idx = phase[i].strand1.start[s]:(phase[i].strand1.start[s + 1] - 1)
            X[idx, i] = H[idx, phase[i].strand1.haplotypelabel[s]]
        end
        idx = phase[i].strand1.start[end]:phase[i].strand1.length
        X[idx, i] = H[idx, phase[i].strand1.haplotypelabel[end]]
        for s in 1:(length(phase[i].strand2.start) - 1)
            idx = phase[i].strand2.start[s]:(phase[i].strand2.start[s + 1] - 1)
            X[idx, i] += H[idx, phase[i].strand2.haplotypelabel[s]]
        end
        idx = phase[i].strand2.start[end]:phase[i].strand2.length
        X[idx, i] += H[idx, phase[i].strand2.haplotypelabel[end]]
    end
end

"""
    impute_discard_phase!(X, H, phase)

Imputes missing entries of `X` using corresponding haplotypes `H` via `phase` information. 
Non-missing entries in `X` will not change, but X and H has to be aligned. This does NOT 
preserve phase information. 
"""
function impute_discard_phase!(
    X::AbstractMatrix,
    compressed_Hunique::CompressedHaplotypes,
    phase::Vector{HaplotypeMosaicPair}
    )

    p, n = size(X)

    @inbounds for person in 1:n, snp in 1:p
        if ismissing(X[snp, person])
            #find which segment the snp is located
            hap1_segment = searchsortedlast(phase[person].strand1.start, snp)
            hap2_segment = searchsortedlast(phase[person].strand2.start, snp)

            #find haplotype pair in corresponding window for this segment
            h1 = phase[person].strand1.haplotypelabel[hap1_segment]
            h2 = phase[person].strand2.haplotypelabel[hap2_segment]
            w1 = phase[person].strand1.window[hap1_segment]
            w2 = phase[person].strand2.window[hap2_segment]
            i1 = snp - compressed_Hunique.start[w1] + 1
            i2 = snp - compressed_Hunique.start[w2] + 1

            # imputation step
            H1 = compressed_Hunique.CW[w1].uniqueH
            H2 = compressed_Hunique.CW[w2].uniqueH
            X[snp, person] = H1[i1, h1] + H2[i2, h2]
        end
    end

    return nothing
end

"""
    update_marker_position!(phaseinfo, tgtfile, reffile)
Converts `phaseinfo`'s strand1 and strand2's starting position in 
terms of matrix rows of `X` to starting position in terms matrix
rows in `H`. 
"""
function update_marker_position!(
    phaseinfo::Vector{HaplotypeMosaicPair},
    XtoH_idx::AbstractVector, 
    )
    people = length(phaseinfo)

    for j in 1:people
        # update strand1's starting position
        for (i, idx) in enumerate(phaseinfo[j].strand1.start)
            phaseinfo[j].strand1.start[i] = XtoH_idx[idx]
        end
        # update strand2's starting position
        for (i, idx) in enumerate(phaseinfo[j].strand2.start)
            phaseinfo[j].strand2.start[i] = XtoH_idx[idx]
        end
    end

    # update first starting position
    for j in 1:people
        phaseinfo[j].strand1.start[1] = 1
        phaseinfo[j].strand2.start[1] = 1
    end

    return nothing
end
