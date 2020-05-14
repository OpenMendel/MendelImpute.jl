"""
    impute_typed_only(tgtfile, reffile, outfile, ph, H, chunks, snps_per_chunk)

Phases and imputes `tgtfile` using `phaseinfo` and outputs result in `outfile`. All genotypes 
in `outfile` are non-missing and phased. 
"""
function impute!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    phaseinfo::Vector{HaplotypeMosaicPair},
    outfile::AbstractString,
    X_sampleID::AbstractVector,
    chr::AbstractVector,
    pos::AbstractVector, 
    ids::AbstractVector,
    ref::AbstractVector,
    alt::AbstractVector
    )
    # impute without changing observed entries
    impute2!(X, H, phaseinfo)

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
    for i in 1:size(X, 1)
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
    update_marker_position!(phaseinfo, tgtfile, reffile)

Converts `phaseinfo`'s strand1 and strand2's starting position in 
terms of matrix rows of `X` to starting position in terms matrix
rows in `H`. 
"""
function update_marker_position!(
    phaseinfo::Vector{HaplotypeMosaicPair},
    XtoH_idx::AbstractVector, 
    ref_records::Int
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
