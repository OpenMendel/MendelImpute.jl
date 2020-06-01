"""
    impute_typed_only(tgtfile, reffile, outfile, ph, H, chunks, snps_per_chunk)

Phases and imputes `tgtfile` using `phaseinfo` and outputs result in `outfile`. All genotypes 
in `outfile` are non-missing and phased. 
"""
function impute!(
    X::AbstractMatrix,
    compressed_haplotypes::CompressedHaplotypes,
    phaseinfo::Vector{HaplotypeMosaicPair},
    outfile::AbstractString,
    X_sampleID::AbstractVector,
    )
    # impute without changing observed entries
    impute_start = time()
    impute2!(X, compressed_haplotypes, phaseinfo)
    impute_time = time() - impute_start

    # retrieve reference file information
    write_start = time()
    chr = compressed_haplotypes.chr
    pos = compressed_haplotypes.pos
    ids = compressed_haplotypes.SNPid
    ref = compressed_haplotypes.refallele
    alt = compressed_haplotypes.altallele

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

    write_time = time() - write_start
    println("impute time  = $(round(impute_time, sigdigits=3)) seconds")
    println("writing time = $(round(write_time, sigdigits=3)) seconds\n")

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
    compressed_haplotypes::CompressedHaplotypes,
    phase::Vector{HaplotypeMosaicPair}
    )

    p, n = size(X)
    width = compressed_haplotypes.width
    CWrange_begin = first.(compressed_haplotypes.CWrange)
    println(typeof(CWrange_begin))
    fdsa

    @inbounds for person in 1:n
        window = 1
        window_counter = 0
        for snp in 1:p
            if ismissing(X[snp, person])
                #find which segment the snp is located
                hap1_segment = searchsortedlast(phase[person].strand1.start, snp)
                hap2_segment = searchsortedlast(phase[person].strand2.start, snp)
                # window = searchsortedlast(CWrange_begin, snp)

                #find haplotype pair in this segment (note: the pair indexes to the entire haplotype pool)
                hap1 = phase[person].strand1.haplotypelabel[hap1_segment]
                hap2 = phase[person].strand2.haplotypelabel[hap2_segment]

                # map hap1 and hap2 back to unique index
                # h1 = complete_idx_to_unique_idx(hap1, window, compressed_haplotypes)
                # h2 = complete_idx_to_unique_idx(hap2, window, compressed_haplotypes)
                h1 = h2 = 1

                # imputation step
                i = snp - (window - 1) * width
                H = compressed_haplotypes[window].uniqueH
                X[snp, person] = H[i, h1] + H[i, h2]
            end
            window_counter += 1
            if window_counter % width == 0 && window != total_windows
                window_counter = 1
            end
        end
    end

    return nothing
end
