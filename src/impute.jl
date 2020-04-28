"""
    impute_typed_only(tgtfile, reffile, outfile, ph, H, chunks, snps_per_chunk)

Phases and imputes `tgtfile` using `phaseinfo` and outputs result in `outfile`. All genotypes 
in `outfile` are non-missing and phased. Markers that are typed in `reffile` but not in 
`tgtfile` will not be in `outfile`. 
"""
function impute_typed_only(
    phaseinfo::Vector{HaplotypeMosaicPair},
    X::AbstractMatrix,
    H_aligned::AbstractMatrix,
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
function impute_untyped(
    tgtfile::AbstractString,
    reffile::AbstractString,
    outfile::AbstractString,
    phaseinfo::Vector{HaplotypeMosaicPair},
    H::AbstractMatrix,
    chunks::Int,
    snps_per_chunk::Int,
    snps_in_last_window::Int
    )

    # some constants
    total_snps = (chunks - 1) * snps_per_chunk + snps_in_last_window
    people = nsamples(tgtfile)
    haplotypes = size(H, 2)
    sample_masks = falses(nsamples(reffile)) # needed for filtering ref_record
    sample_masks[1:people] .= true

    # convert phase's starting position from matrix index to marker position
    update_marker_position!(phaseinfo, tgtfile)

    # write phase information to outfile
    tgt_reader = VCF.Reader(openvcf(tgtfile, "r"))
    ref_reader = VCF.Reader(openvcf(reffile, "r"))
    tgt_record = read(tgt_reader) # first record
    tgt_pos = VCF.pos(tgt_record) # first record's position
    writer = VCF.Writer(openvcf(outfile, "w"), header(tgt_reader))

    if chunks > 1
        # reassign and update H chunk by chunk
        # Hreader = VCF.Reader(openvcf(reffile, "r"))
        # H = BitArray{2}(undef, snps_per_chunk, haplotypes)
        # copy_ht_trans!(H, Hreader)
        # record_counter = chunk_counter = 1
        # pmeter = Progress(total_snps, 5, "Writing to file...")

        # for (i, ref_record) in enumerate(ref_reader)
        #     ref_pos = VCF.pos(ref_record)
        #     if ref_pos != tgt_pos
        #         gtkey = VCF.findgenokey(ref_record, "GT")
        #         if !isnothing(gtkey) 
        #             # filter record so it only contains as many people as in target
        #             VCFTools.filter_record!(ref_record, sample_masks)

        #             # if snp exist only in reference file, fetch nearest haplotypelabel 
        #             for (person, geno) in enumerate(ref_record.genotype)
        #                 #find where snp is located in phase. max() avoids indexing error when ref snps occurs before 1st target snp
        #                 hap1_position = max(1, searchsortedlast(phaseinfo[person].strand1.start, ref_pos))
        #                 hap2_position = max(1, searchsortedlast(phaseinfo[person].strand2.start, ref_pos))

        #                 #find the correct haplotypes 
        #                 hap1 = phaseinfo[person].strand1.haplotypelabel[hap1_position]
        #                 hap2 = phaseinfo[person].strand2.haplotypelabel[hap2_position]

        #                 # save actual allele to data. "0" (REF) => 0x30, "1" (ALT) => 0x31
        #                 row = i - (chunk_counter - 1) * snps_per_chunk
        #                 a1, a2 = H[row, hap1], H[row, hap2]
        #                 ref_record.data[geno[gtkey][1]] = ifelse(a1, 0x31, 0x30)
        #                 ref_record.data[geno[gtkey][2]] = 0x7c # phased data has separator '|'
        #                 ref_record.data[geno[gtkey][3]] = ifelse(a2, 0x31, 0x30)
        #             end
        #             write(writer, ref_record)
        #         end
        #     else
        #         gtkey = VCF.findgenokey(tgt_record, "GT")
        #         if !isnothing(gtkey) 
        #             # if snp exist in target, loop over samples and change only missing entries
        #             for (person, geno) in enumerate(tgt_record.genotype)
        #                 if tgt_record.data[geno[gtkey][1]] == 0x2e # 0x2e is '.' which indicates missing
        #                     #find where snp is located in phase
        #                     hap1_position = searchsortedlast(phaseinfo[person].strand1.start, tgt_pos)
        #                     hap2_position = searchsortedlast(phaseinfo[person].strand2.start, tgt_pos)

        #                     #find the correct haplotypes 
        #                     hap1 = phaseinfo[person].strand1.haplotypelabel[hap1_position]
        #                     hap2 = phaseinfo[person].strand2.haplotypelabel[hap2_position]

        #                     # save actual allele to data. "0" (REF) => 0x30, "1" (ALT) => 0x31
        #                     row = i - (chunk_counter - 1) * snps_per_chunk
        #                     a1, a2 = H[row, hap1], H[row, hap2]
        #                     tgt_record.data[geno[gtkey][1]] = ifelse(a1, 0x31, 0x30)
        #                     tgt_record.data[geno[gtkey][2]] = 0x7c # phased data has separator '|'
        #                     tgt_record.data[geno[gtkey][3]] = ifelse(a2, 0x31, 0x30)
        #                 end
        #             end
        #             write(writer, tgt_record)
        #         end

        #         # read next target record
        #         if !eof(tgt_reader)
        #             tgt_record = read(tgt_reader) 
        #             tgt_pos = VCF.pos(tgt_record)
        #         end
        #     end

        #     # move to next chunk if we reached the end of current chunk 
        #     record_counter += 1
        #     if record_counter > snps_per_chunk
        #         chunk_counter += 1
        #         record_counter = 1
        #         chunk_counter == chunks && (H = BitArray{2}(undef, snps_in_last_window, haplotypes)) #resize H
        #         copy_ht_trans!(H, Hreader)
        #     end

        #     # update progress
        #     next!(pmeter) 
        # end
        # close(Hreader)
    else
        # loop over each record (snp) in ref file
        pmeter = Progress(size(H, 1), 5, "Writing to file...")
        for (i, ref_record) in enumerate(ref_reader)
            ref_pos = VCF.pos(ref_record)
            if ref_pos != tgt_pos
                # if snp exist only in reference file, fetch nearest haplotypelabel for everybody
                gtkey = VCF.findgenokey(ref_record, "GT")
                if !isnothing(gtkey) 
                    # filter record so it only contains as many people as in target
                    VCFTools.filter_record!(ref_record, sample_masks)

                    for (person, geno) in enumerate(ref_record.genotype)
                        #find where snp is located in phase. max() avoids indexing error when ref snps occurs before 1st target snp
                        hap1_position = max(1, searchsortedlast(phaseinfo[person].strand1.start, ref_pos))
                        hap2_position = max(1, searchsortedlast(phaseinfo[person].strand2.start, ref_pos))

                        #find the correct haplotypes 
                        hap1 = phaseinfo[person].strand1.haplotypelabel[hap1_position]
                        hap2 = phaseinfo[person].strand2.haplotypelabel[hap2_position]

                        # save actual allele to data. "0" (REF) => 0x30, "1" (ALT) => 0x31
                        a1, a2 = H[i, hap1], H[i, hap2]
                        ref_record.data[geno[gtkey][1]] = ifelse(a1, 0x31, 0x30)
                        ref_record.data[geno[gtkey][2]] = 0x7c # phased data has separator '|'
                        ref_record.data[geno[gtkey][3]] = ifelse(a2, 0x31, 0x30)
                    end
                    write(writer, ref_record)
                end
            else
                # if snp exist in target, loop over samples and change only missing entries
                gtkey = VCF.findgenokey(tgt_record, "GT")
                if !isnothing(gtkey) 
                    for (person, geno) in enumerate(tgt_record.genotype)
                        if tgt_record.data[geno[gtkey][1]] == 0x2e # 0x2e is '.' which indicates missing
                            #find where snp is located in phase
                            hap1_position = searchsortedlast(phaseinfo[person].strand1.start, tgt_pos)
                            hap2_position = searchsortedlast(phaseinfo[person].strand2.start, tgt_pos)

                            #find the correct haplotypes 
                            hap1 = phaseinfo[person].strand1.haplotypelabel[hap1_position]
                            hap2 = phaseinfo[person].strand2.haplotypelabel[hap2_position]

                            # save actual allele to data. "0" (REF) => 0x30, "1" (ALT) => 0x31
                            a1, a2 = H[i, hap1], H[i, hap2]
                            tgt_record.data[geno[gtkey][1]] = ifelse(a1, 0x31, 0x30)
                            tgt_record.data[geno[gtkey][2]] = 0x7c # phased data has separator '|'
                            tgt_record.data[geno[gtkey][3]] = ifelse(a2, 0x31, 0x30)
                        end
                    end
                    write(writer, tgt_record)
                end

                # read next target record
                if !eof(tgt_reader)
                    tgt_record = read(tgt_reader) 
                    tgt_pos = VCF.pos(tgt_record)
                end
            end

            next!(pmeter) #update progress
        end
    end

    # close reader/writers 
    flush(writer); close(writer); close(tgt_reader); close(ref_reader)
end

function impute_untyped2(
    tgtfile::AbstractString,
    reffile::AbstractString,
    outfile::AbstractString,
    phaseinfo::Vector{HaplotypeMosaicPair},
    H::AbstractMatrix,
    chunks::Int,
    snps_per_chunk::Int,
    snps_in_last_window::Int
    )
    
    # some constants and readers
    people = nsamples(tgtfile)
    ref_reader = VCF.Reader(openvcf(reffile, "r"))
    tgt_reader = VCF.Reader(openvcf(tgtfile, "r"))

    # convert phase's starting position from tgt matrix index to ref matrix index
    update_marker_position!(phaseinfo, tgtfile, reffile)

    if chunks > 1
        # TODO
    else
        # first impute
        X = zeros(Int, size(H, 1), people)
        impute!(X, H, phaseinfo)

        # write minimal meta information to outfile
        io = openvcf(outfile, "w")
        write(io, "##fileformat=VCFv4.2\n")
        write(io, "##source=MendelImpute\n")
        write(io, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

        # header line should match reffile (i.e. sample ID's should match)
        write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for id in VCF.header(tgt_reader).sampleID
            write(io, "\t", id)
        end
        write(io, "\n")
        close(tgt_reader)

        # write phase info
        pmeter = Progress(size(X, 1), 5, "Writing to file...")
        for snp in 1:size(X, 1)
            ref_record = read(ref_reader)
            write(io, VCF.chrom(ref_record), "\t", string(VCF.pos(ref_record)), "\t", 
                VCF.id(ref_record)[1], "\t", VCF.ref(ref_record), "\t", VCF.alt(ref_record)[1], 
                "\t.\tPASS\t.\tGT")
            for i in 1:size(X, 2)
                if X[snp, i] == 2
                    write(io, "\t1/1")
                elseif X[snp, i] == 1
                    write(io, "\t1/0")
                else
                    write(io, "\t0/0")
                end
            end
            write(io, "\n")
            next!(pmeter)
        end
        close(io); close(ref_reader)
    end
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
terms of matrix rows of `tgtfile` to starting position in terms matrix
rows in `reffile`. 

TODO: iterating over reference vcf files takes a long time, can we get POS somewhere else?
"""
function update_marker_position!(
    phaseinfo::Vector{HaplotypeMosaicPair},
    tgtfile::AbstractString,
    reffile::AbstractString
    )
    people = length(phaseinfo)
    tgt_reader = VCF.Reader(openvcf(tgtfile, "r"))
    ref_reader = VCF.Reader(openvcf(reffile, "r"))
    ref_records = nrecords(reffile)

    # find marker position for each SNP in reference file
    # TODO: this loop is extremely slow
    ref_marker_pos = zeros(Int, ref_records)
    for (i, record) in enumerate(ref_reader)
        ref_marker_pos[i] = VCF.pos(record)
    end
    close(ref_reader)

    # find marker position for each SNP in target file
    tgt_marker_pos = zeros(Int, phaseinfo[1].strand1.length)
    for (i, record) in enumerate(tgt_reader)
        tgt_marker_pos[i] = VCF.pos(record)
    end
    close(tgt_reader)

    for j in 1:people
        # update strand1's starting position (optimization: can change findfirst to findnext)
        for (i, idx) in enumerate(phaseinfo[j].strand1.start)
            phaseinfo[j].strand1.start[i] = findfirst(x -> x == tgt_marker_pos[idx], ref_marker_pos) :: Int
        end
        # update strand2's starting position
        for (i, idx) in enumerate(phaseinfo[j].strand2.start)
            phaseinfo[j].strand2.start[i] = findfirst(x -> x == tgt_marker_pos[idx], ref_marker_pos) :: Int
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
