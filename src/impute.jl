"""
    impute_typed_only(tgtfile, reffile, outfile, ph, H, chunks, snps_per_chunk)

Phases and imputes `tgtfile` using `phaseinfo` and outputs result in `outfile`. All genotypes 
in `outfile` are non-missing and phased. Markers that are typed in `reffile` but not in 
`tgtfile` will not be in `outfile`. 
"""
function impute_typed_only!(
    X::AbstractMatrix,
    H_aligned::AbstractMatrix,
    phaseinfo::Vector{HaplotypeMosaicPair},
    outfile::AbstractString,
    X_sampleID::AbstractVector,
    X_chr::AbstractVector,
    X_pos::AbstractVector, 
    X_ids::AbstractVector,
    X_ref::AbstractVector,
    X_alt::AbstractVector
    )
    # impute without changing observed entries
    impute2!(X, H_aligned, phaseinfo)

    # write minimal meta information to outfile
    io = openvcf(outfile, "w")
    write(io, "##fileformat=VCFv4.2\n")
    write(io, "##source=MendelImpute\n")
    write(io, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

    # header line should match reffile (i.e. sample ID's should match)
    write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    for id in X_sampleID
        write(io, "\t", id)
    end
    write(io, "\n")

    pmeter = Progress(size(X, 1), 5, "Writing to file...")
    for i in 1:size(X, 1)
        # write meta info (chrom/pos/id/ref/alt)
        write(io, X_chr[i], "\t", string(X_pos[i]), "\t", X_ids[i][1], "\t", X_ref[i], "\t", X_alt[i][1], "\t.\tPASS\t.\tGT")
        for j in 1:size(X, 2)
            if X[i, j] == 0
                write(io, "\t0/0")
            elseif X[i, j] == 1
                write(io, "\t1/0")
            else
                write(io, "\t1/1")
            end
        end
        write(io, "\n")
        next!(pmeter)
    end

    # close & return
    close(io)
    return X
end

"""
    impute_untyped(tgtfile, reffile, outfile, ph, H, chunks, snps_per_chunk)

Phases and imputes `tgtfile` using `phaseinfo` and outputs result in `outfile`. All genotypes 
in `outfile` are non-missing and phased. Markers that are typed in `reffile` but not in 
`tgtfile` (determined via SNP position) will be imputed in `outfile` as well. 
"""
function impute_untyped!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    phaseinfo::Vector{HaplotypeMosaicPair},
    outfile::AbstractString,
    X_sampleID::AbstractVector,
    H_pos::AbstractVector, 
    H_chr::AbstractVector,
    H_ids::AbstractVector,
    H_ref::AbstractVector,
    H_alt::AbstractVector
    )
    # impute without changing observed entries
    impute2!(X, H, phaseinfo)

    # write minimal meta information to outfile
    io = openvcf(outfile, "w")
    write(io, "##fileformat=VCFv4.2\n")
    write(io, "##source=MendelImpute\n")
    write(io, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

    # header line should match reffile (i.e. sample ID's should match)
    write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    for id in X_sampleID
        write(io, "\t", id)
    end
    write(io, "\n")

    pmeter = Progress(size(H, 1), 5, "Writing to file...")
    for i in 1:size(X, 1)
        # write meta info (chrom/pos/id/ref/alt)
        write(io, H_chr[i], "\t", string(H_pos[i]), "\t", H_ids[i][1], "\t", H_ref[i], "\t", H_alt[i][1], "\t.\tPASS\t.\tGT")
        for j in 1:size(X, 2)
            if X[i, j] == 0
                write(io, "\t0/0")
            elseif X[i, j] == 1
                write(io, "\t1/0")
            else
                write(io, "\t1/1")
            end
        end
        write(io, "\n")
        next!(pmeter)
    end

    # close & return
    close(io)
    return X
end

"""
    update_marker_position!(phaseinfo, tgtfile)

Converts `phaseinfo`'s strand1 and strand2's starting position in 
terms of matrix rows to starting position in terms of SNP position. 

TODO: iterating over vcf files takes a long time, can we get POS somewhere else?
"""
function update_marker_position!(
    phaseinfo::Vector{HaplotypeMosaicPair},
    tgtfile::AbstractString;
    )
    people = length(phaseinfo)
    reader = VCF.Reader(openvcf(tgtfile, "r"))
    marker_pos = zeros(Int, phaseinfo[1].strand1.length)

    # find marker position for each SNP
    for (i, record) in enumerate(reader)
        gtkey = VCF.findgenokey(record, "GT")
        if !isnothing(gtkey) 
            marker_pos[i] = VCF.pos(record)
        end
    end

    for j in 1:people
        # update strand1's starting position
        for (i, idx) in enumerate(phaseinfo[j].strand1.start)
            phaseinfo[j].strand1.start[i] = marker_pos[idx]
        end
        # update strand2's starting position
        for (i, idx) in enumerate(phaseinfo[j].strand2.start)
            phaseinfo[j].strand2.start[i] = marker_pos[idx]
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
    tgt_marker_pos::AbstractVector, 
    ref_marker_pos::AbstractVector
    )
    people = length(phaseinfo)
    ref_records = length(ref_marker_pos)

    # for j in 1:people
    #     # update strand1's starting position (optimization: can change findfirst to findnext)
    #     for (i, idx) in enumerate(phaseinfo[j].strand1.start)
    #         phaseinfo[j].strand1.start[i] = findfirst(x -> x == tgt_marker_pos[idx], ref_marker_pos) :: Int
    #     end
    #     # update strand2's starting position
    #     for (i, idx) in enumerate(phaseinfo[j].strand2.start)
    #         phaseinfo[j].strand2.start[i] = findfirst(x -> x == tgt_marker_pos[idx], ref_marker_pos) :: Int
    #     end
    # end

    for j in 1:people
        # update strand1's starting position
        for i in eachindex(phaseinfo[j].strand1.start)
            phaseinfo[j].strand1.start[i] = XtoH_idx[i]
        end
        # update strand2's starting position
        for i in eachindex(phaseinfo[j].strand2.start)
            phaseinfo[j].strand2.start[i] = XtoH_idx[i]
        end
    end

    # update first starting position and length
    for j in 1:people
        phaseinfo[j].strand1.start[1] = 1
        phaseinfo[j].strand2.start[1] = 1
        phaseinfo[j].strand1.length = ref_records
        phaseinfo[j].strand2.length = ref_records
    end

    return nothing
end

"""
    impute!(X, H, phase)

Imputes `X` completely using segments of haplotypes `H` where segments are stored in `phase`. 
Non-missing entries in `X` can be different after imputation. 
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
    impute2!(X, H, phase)

Imputes missing entries of `X` using corresponding haplotypes `H` via `phase` information. 
Non-missing entries in `X` will not change, but X and H has to be aligned. 
"""
function impute2!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    phase::Vector{HaplotypeMosaicPair}
    )

    p, n = size(X)

    @inbounds for snp in 1:p, person in 1:n
        if ismissing(X[snp, person])
            #find where snp is located in phase
            hap1_position = searchsortedlast(phase[person].strand1.start, snp)
            hap2_position = searchsortedlast(phase[person].strand2.start, snp)

            #find the correct haplotypes 
            hap1 = phase[person].strand1.haplotypelabel[hap1_position]
            hap2 = phase[person].strand2.haplotypelabel[hap2_position]

            # imputation step 
            X[snp, person] = H[snp, hap1] + H[snp, hap2]
        end
    end

    return nothing
end
