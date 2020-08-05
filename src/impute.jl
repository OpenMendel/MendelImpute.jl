"""
    output_unphased!(X, compressed_haplotypes, phaseinfo, outfile, X_sampleID)

Imputes `X` using `phaseinfo` and outputs result in `outfile`. All genotypes
in `outfile` are non-missing and unphased.

If `XtoH_idx == nothing`, all SNPs in reference file will be imputed.
"""
function Base.write(
    outfile::AbstractString,
    X::AbstractMatrix,
    compressed_haplotypes::CompressedHaplotypes,
    X_sampleID::AbstractVector,
    XtoH_idx::Union{Nothing, AbstractVector} = nothing,
    )
    # retrieve reference file information
    chr = (isnothing(XtoH_idx) ? compressed_haplotypes.chr : 
                                 compressed_haplotypes.chr[XtoH_idx])
    pos = (isnothing(XtoH_idx) ? compressed_haplotypes.pos : 
                                 compressed_haplotypes.pos[XtoH_idx])
    ids = (isnothing(XtoH_idx) ? compressed_haplotypes.SNPid : 
                                 compressed_haplotypes.SNPid[XtoH_idx])
    ref = (isnothing(XtoH_idx) ? compressed_haplotypes.refallele : 
                                 compressed_haplotypes.refallele[XtoH_idx])
    alt = (isnothing(XtoH_idx) ? compressed_haplotypes.altallele : 
                                 compressed_haplotypes.altallele[XtoH_idx])

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
        # write meta info (chrom/pos/snpid/ref/alt)
        print(pb, chr[i], "\t", string(pos[i]), "\t", ids[i][1], "\t", ref[i],
            "\t", alt[i][1], "\t.\tPASS\t.\tGT")

        # print ith record
        write_snp!(pb, @view(X[i, :]))

        (bytesavailable(pb) > (16*1024)) && write(io, take!(pb))
        next!(pmeter)
    end
    write(io, take!(pb))

    # close & return
    close(io); close(pb)
    return nothing
end

"""
    write(outfile, X1, X2, compressed_haplotypes, phaseinfo, X_sampleID)

All genotypes in `outfile` are non-missing and phased.
"""
function Base.write(
    outfile::AbstractString,
    X1::AbstractMatrix,
    X2::AbstractMatrix,
    compressed_haplotypes::CompressedHaplotypes,
    X_sampleID::AbstractVector;
    XtoH_idx::Union{Nothing, AbstractVector} = nothing,
    )
    # retrieve reference file information
    chr = (isnothing(XtoH_idx) ? compressed_haplotypes.chr :
                                 compressed_haplotypes.chr[XtoH_idx])
    pos = (isnothing(XtoH_idx) ? compressed_haplotypes.pos :
                                 compressed_haplotypes.pos[XtoH_idx])
    ids = (isnothing(XtoH_idx) ? compressed_haplotypes.SNPid :
                                 compressed_haplotypes.SNPid[XtoH_idx])
    ref = (isnothing(XtoH_idx) ? compressed_haplotypes.refallele :
                                 compressed_haplotypes.refallele[XtoH_idx])
    alt = (isnothing(XtoH_idx) ? compressed_haplotypes.altallele :
                                 compressed_haplotypes.altallele[XtoH_idx])

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

    pmeter = Progress(size(X1, 1), 5, "Writing to file...")
    @inbounds for i in 1:size(X1, 1)
        # write meta info (chrom/pos/id/ref/alt)
        print(pb, chr[i], "\t", string(pos[i]), "\t", ids[i][1], "\t", ref[i],
            "\t", alt[i][1], "\t.\tPASS\t.\tGT")

        # print ith record
        write_snp!(pb, @view(X1[i, :]), @view(X2[i, :]))

        (bytesavailable(pb) > (16*1024)) && write(io, take!(pb))
        next!(pmeter)
    end
    write(io, take!(pb))

    # close & return
    close(io); close(pb)
    return nothing
end

"""
Helper function for saving a record (SNP), not tracking phase information.
"""
function write_snp!(pb::IOBuffer, X::AbstractVector)
    n = length(X)
    @inbounds for j in 1:n
        if X[j] == 0
            print(pb, "\t0/0")
        elseif X[j] == 1
            print(pb, "\t1/0")
        elseif X[j] == 2
            print(pb, "\t1/1")
        else
            error("imputed genotypes can only be 0, 1, 2 but got $(X[j])")
        end
    end
    print(pb, "\n")
    nothing
end

"""
Helper function for saving a record (SNP), tracking phase information.
"""
function write_snp!(pb::IOBuffer, X1::AbstractVector, X2::AbstractVector)
    n = length(X1)
    @assert n == length(X2)
    @inbounds for j in 1:n
        if X1[j] == X2[j] == 0
            print(pb, "\t0|0")
        elseif X1[j] == 0 && X2[j] == 1
            print(pb, "\t0|1")
        elseif X1[j] == 1 && X2[j] == 0
            print(pb, "\t1|0")
        elseif X1[j] == 1 && X2[j] == 1
            print(pb, "\t1|1")
        else
            error("phased genotypes can only be 0|0, 0|1, 1|0 or 1|1 but
                got $(X1[j])|$(X2[j])")
        end
    end
    print(pb, "\n")
    nothing
end

"""
    impute!(X1, X2, H, phase)

Imputes `X = X1 + X2` completely using haplotype segments of `H`, where segments
information are stored in `phase`. `X1` is strand1 and `X2` is strand 2.

Non-missing entries in `X` can be different after imputation.
"""
function impute!(
    X1::AbstractMatrix,
    X2::AbstractMatrix,
    compressed_Hunique::CompressedHaplotypes,
    phase::Vector{HaplotypeMosaicPair}
    )

    fill!(X1, 0)
    fill!(X2, 0)

    # loop over individuals
    for i in 1:size(X1, 2)
        # strand 1
        for s in 1:(length(phase[i].strand1.start) - 1)
            X_idx = phase[i].strand1.start[s]:(phase[i].strand1.start[s + 1] - 1)
            w = phase[i].strand1.window[s]
            H = compressed_Hunique.CW[w].uniqueH
            H_start = abs(phase[i].strand1.start[s] - 
                compressed_Hunique.Hstart[w]) + 1
            H_idx = H_start:(H_start + length(X_idx) - 1)
            X1[X_idx, i] = H[H_idx, phase[i].strand1.haplotypelabel[s]]
        end
        w = phase[i].strand1.window[end]
        X_idx = phase[i].strand1.start[end]:phase[i].strand1.length
        H_start = abs(phase[i].strand1.start[end] - 
            compressed_Hunique.Hstart[w]) + 1
        H_idx = H_start:(H_start + length(X_idx) - 1)
        H = compressed_Hunique.CW[w].uniqueH
        X1[X_idx, i] = H[H_idx, phase[i].strand1.haplotypelabel[end]]

        # strand 2
        for s in 1:(length(phase[i].strand2.start) - 1)
            X_idx = phase[i].strand2.start[s]:(phase[i].strand2.start[s + 1] - 1)
            w = phase[i].strand2.window[s]
            H = compressed_Hunique.CW[w].uniqueH
            H_start = abs(phase[i].strand2.start[s] - 
                compressed_Hunique.Hstart[w]) + 1
            H_idx = H_start:(H_start + length(X_idx) - 1)
            X2[X_idx, i] = H[H_idx, phase[i].strand2.haplotypelabel[s]]
        end
        X_idx = phase[i].strand2.start[end]:phase[i].strand2.length
        w = phase[i].strand2.window[end]
        H = compressed_Hunique.CW[w].uniqueH
        H_start = abs(phase[i].strand2.start[end] - 
            compressed_Hunique.Hstart[w]) + 1
        H_idx = H_start:(H_start + length(X_idx) - 1)
        X2[X_idx, i] = H[H_idx, phase[i].strand2.haplotypelabel[end]]
    end
end

"""
    impute_discard_phase!(X, H, phase)

Imputes missing entries of `X` using corresponding haplotypes `H` via `phase`
information. Non-missing entries in `X` will not change, but X and H has to be
aligned. This does NOT preserve phase information.
"""
function impute_discard_phase!(
    X::AbstractMatrix,
    compressed_Hunique::CompressedHaplotypes,
    phase::Vector{HaplotypeMosaicPair}
    )

    p, n = size(X)

    # for person in 1:n
    ThreadPools.@qthreads for person in 1:n
        @inbounds for snp in 1:p
            if ismissing(X[snp, person])
                #find which segment the snp is located
                hap1_segment = searchsortedlast(phase[person].strand1.start,snp)
                hap2_segment = searchsortedlast(phase[person].strand2.start,snp)

                #find haplotype pair in corresponding window for this segment
                h1 = phase[person].strand1.haplotypelabel[hap1_segment]
                h2 = phase[person].strand2.haplotypelabel[hap2_segment]
                w1 = phase[person].strand1.window[hap1_segment]
                w2 = phase[person].strand2.window[hap2_segment]
                i1 = snp - compressed_Hunique.Hstart[w1] + 1
                i2 = snp - compressed_Hunique.Hstart[w2] + 1

                # imputation step
                try
                    H1 = compressed_Hunique.CW[w1].uniqueH
                    H2 = compressed_Hunique.CW[w2].uniqueH
                    X[snp, person] = H1[i1, h1] + H2[i2, h2]
                catch
                    println("snp = $snp, h1 = $h1, w1 = $w1, i1 = $i1")
                    println("snp = $snp, h2 = $h2, w2 = $w2, i2 = $i2\n")

                    break
                end
            end
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

    @inbounds for j in 1:people
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
    @inbounds for j in 1:people
        phaseinfo[j].strand1.start[1] = 1
        phaseinfo[j].strand2.start[1] = 1
    end

    return nothing
end
