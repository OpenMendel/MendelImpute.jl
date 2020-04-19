"""
    impute_typed_only(tgtfile, reffile, outfile, ph, H, chunks, snps_per_chunk)

Phases and imputes `tgtfile` using `phaseinfo` and outputs result in `outfile`. All genotypes 
in `outfile` are non-missing and phased. Markers that are typed in `reffile` but not in 
`tgtfile` will not be in `outfile`. 
"""
function impute_typed_only(
    tgtfile::AbstractString,
    reffile::AbstractString,
    outfile::AbstractString,
    phaseinfo::Vector{HaplotypeMosaicPair},
    H::AbstractMatrix,
    chunks::Int,
    snps_per_chunk::Int,
    )
    haplotypes = size(H, 2)

    # write phase information to outfile
    reader = VCF.Reader(openvcf(tgtfile, "r"))
    writer = VCF.Writer(openvcf(outfile, "w"), header(reader))
    pmeter = Progress(nrecords(tgtfile), 5, "Writing to file...")
    if chunks > 1
        # reassign and update H chunk by chunk
        Hreader = VCF.Reader(openvcf(reffile, "r"))
        H = BitArray{2}(undef, snps_per_chunk, haplotypes)
        copy_ht_trans!(H, Hreader)
        record_counter = chunk_counter = 1
        for (i, record) in enumerate(reader)
            gtkey = VCF.findgenokey(record, "GT")
            if !isnothing(gtkey) 
                # loop over samples
                for (j, geno) in enumerate(record.genotype)
                    # if missing = '.' = 0x2e
                    if record.data[geno[gtkey][1]] == 0x2e
                        #find where snp is located in phase
                        hap1_position = searchsortedlast(phaseinfo[j].strand1.start, i)
                        hap2_position = searchsortedlast(phaseinfo[j].strand2.start, i)

                        #find the correct haplotypes 
                        hap1 = phaseinfo[j].strand1.haplotypelabel[hap1_position]
                        hap2 = phaseinfo[j].strand2.haplotypelabel[hap2_position]

                        # save actual allele to data. "0" (REF) => 0x30, "1" (ALT) => 0x31
                        row = i - (chunk_counter - 1) * snps_per_chunk
                        a1, a2 = H[row, hap1], H[row, hap2]
                        record.data[geno[gtkey][1]] = ifelse(a1, 0x31, 0x30)
                        record.data[geno[gtkey][2]] = 0x7c # phased data has separator '|'
                        record.data[geno[gtkey][3]] = ifelse(a2, 0x31, 0x30)
                    end
                end
            end

            write(writer, record)

            # move to next chunk if we reached the end of current chunk 
            record_counter += 1
            if record_counter > snps_per_chunk
                chunk_counter += 1
                record_counter = 1
                chunk_counter == chunks && (H = BitArray{2}(undef, remaining_snps, haplotypes)) #resize H
                copy_ht_trans!(H, Hreader)
            end

            # update progress
            next!(pmeter) 
        end
        close(Hreader)
    else
        # loop over each record (snp)
        for (i, record) in enumerate(reader)
            gtkey = VCF.findgenokey(record, "GT")
            if !isnothing(gtkey) 
                # loop over samples
                for (j, geno) in enumerate(record.genotype)
                    # if missing = '.' = 0x2e
                    if record.data[geno[gtkey][1]] == 0x2e
                        #find where snp is located in phase
                        hap1_position = searchsortedlast(phaseinfo[j].strand1.start, i)
                        hap2_position = searchsortedlast(phaseinfo[j].strand2.start, i)

                        #find the correct haplotypes 
                        hap1 = phaseinfo[j].strand1.haplotypelabel[hap1_position]
                        hap2 = phaseinfo[j].strand2.haplotypelabel[hap2_position]

                        # save actual allele to data. "0" (REF) => 0x30, "1" (ALT) => 0x31
                        a1, a2 = H[i, hap1], H[i, hap2]
                        record.data[geno[gtkey][1]] = ifelse(a1, 0x31, 0x30)
                        record.data[geno[gtkey][2]] = 0x7c # phased data has separator '|'
                        record.data[geno[gtkey][3]] = ifelse(a2, 0x31, 0x30)
                    end
                end
            end
            write(writer, record)
            next!(pmeter) #update progress
        end
    end

    # close 
    flush(writer); close(reader); close(writer)
end

"""
    impute_untyped(tgtfile, reffile, outfile, ph, H, chunks, snps_per_chunk)

Phases and imputes `tgtfile` using `phaseinfo` and outputs result in `outfile`. All genotypes 
in `outfile` are non-missing and phased. Markers that are typed in `reffile` but not in 
`tgtfile` will be imputed in `outfile` as well. 
"""
function impute_untyped(
    tgtfile::AbstractString,
    reffile::AbstractString,
    outfile::AbstractString,
    phaseinfo::Vector{HaplotypeMosaicPair},
    H::AbstractMatrix,
    chunks::Int,
    snps_per_chunk::Int,
    haplotypes::Int,
    )
    # write phase information to outfile
    reader = VCF.Reader(openvcf(tgtfile, "r"))
    writer = VCF.Writer(openvcf(outfile, "w"), header(reader))
    pmeter = Progress(nrecords(tgtfile), 5, "Writing to file...")
    if chunks > 1
        # TODO
    else
        # loop over each record (snp)
        for record in reader
            pos   = VCF.pos(record)
            gtkey = VCF.findgenokey(record, "GT")
            if !isnothing(gtkey) 
                # loop over samples
                for (j, geno) in enumerate(record.genotype)
                    # if missing = '.' = 0x2e
                    if record.data[geno[gtkey][1]] == 0x2e
                        #find where snp is located in phase
                        hap1_position = searchsortedlast(phaseinfo[j].strand1.start, i)
                        hap2_position = searchsortedlast(phaseinfo[j].strand2.start, i)

                        #find the correct haplotypes 
                        hap1 = phaseinfo[j].strand1.haplotypelabel[hap1_position]
                        hap2 = phaseinfo[j].strand2.haplotypelabel[hap2_position]

                        # save actual allele to data. "0" (REF) => 0x30, "1" (ALT) => 0x31
                        a1, a2 = H[i, hap1], H[i, hap2]
                        record.data[geno[gtkey][1]] = ifelse(a1, 0x31, 0x30)
                        record.data[geno[gtkey][2]] = 0x7c # phased data has separator '|'
                        record.data[geno[gtkey][3]] = ifelse(a2, 0x31, 0x30)
                    end
                end
            end
            write(writer, record)
            next!(pmeter) #update progress
        end
    end

    # close 
    flush(writer); close(reader); close(writer)
end

function update_marker_position!(
    phaseinfo::Vector{HaplotypeMosaicPair},
    tgtfile::AbstractString,
    )
    people = length(phaseinfo)
    reader = VCF.Reader(openvcf(tgtfile, "r"))
    s1_counters, s2_counters = ones(Int, people), ones(Int, people)
    for (i, record) in enumerate(reader)
        pos   = VCF.pos(record)
        gtkey = VCF.findgenokey(record, "GT")
        @inbounds for j in 1:people
            if phaseinfo[j].strand1.start[s1_counters[j]] == i && !isnothing(gtkey) 
                phaseinfo[j].strand1.start[s1_counters[j]] = pos
                s1_counters[j] += 1
            end
            if phaseinfo[j].strand2.start[s2_counters[j]] == i && !isnothing(gtkey) 
                phaseinfo[j].strand2.start[s2_counters[j]] = pos
                s2_counters[j] += 1
            end
        end
    end
end