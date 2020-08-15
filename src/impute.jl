"""
    write(outfile, X, compressed_haplotypes, X_sampleID, XtoH_idx)

Writes imputed `X` into `outfile`. All genotypes in `outfile` are non-missing.
and unphased. 

# Notes
Here the writing routine is emulating `write_dlm` in Base at 
https://github.com/JuliaLang/julia/blob/3608c84e6093594fe86923339fc315231492484c/stdlib/DelimitedFiles/src/DelimitedFiles.jl#L736
"""
function Base.write(
    outfile::AbstractString,
    X::Union{AbstractMatrix, Tuple{AbstractMatrix, AbstractMatrix}},
    compressed_haplotypes::CompressedHaplotypes,
    X_sampleID::AbstractVector,
    XtoH_idx::Union{Nothing, AbstractVector} = nothing,
    )
    threads = Threads.nthreads()
    snps = typeof(X) <: AbstractMatrix ? size(X, 1) : size(X[1], 1)
    len = div(snps, threads)
    files = ["tmp$i.vcf.gz" for i in 1:threads]

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
    io = [openvcf(files[i], "w") for i in 1:threads]
    pb = [PipeBuffer() for _ in 1:threads]
    print(pb[1], "##fileformat=VCFv4.2\n")
    print(pb[1], "##source=MendelImpute\n")
    print(pb[1], "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

    # header line should match reffile (i.e. sample ID's should match)
    print(pb[1], "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    for id in X_sampleID
        print(pb[1], "\t", id)
    end
    print(pb[1], "\n")
    bytesavailable(pb[1]) > 1048576 && write(io[1], take!(pb[1]))

    # each thread writes `len` SNPs
    pmeter = Progress(snps, 5, "Writing to file...")
    Threads.@threads for t in 1:threads
        id = Threads.threadid()
        cur_ranges = (id == threads ? 
            ((threads-1)*len+1:snps) : (1:len) .+ (t-1)*len)

        @inbounds for i in cur_ranges
            # write meta info (chrom/pos/snpid/ref/alt)
            print(pb[id], chr[i], "\t", string(pos[i]), "\t", ids[i][1], "\t", 
                ref[i], "\t", alt[i][1], "\t.\tPASS\t.\tGT")
            # print ith record
            write_snp!(pb[id], X, i) 
            bytesavailable(pb[id]) > 1048576 && write(io[id], take!(pb[id]))
            next!(pmeter)
        end
        write(io[id], take!(pb[id]))
    end
    close.(io); close.(pb) # close io and buffer

    # concatenate all files into 1 VCF file
    run(pipeline(`cat $files`, stdout=outfile))

    # delete intermediate files
    for i in 1:threads
        rm("tmp$i.vcf.gz", force=true)
    end

    return nothing
end

"""
Helper function for saving a record (SNP), not tracking phase information.
"""
function write_snp!(pb::IOBuffer, X::AbstractMatrix, i::Int)
    x = @view(X[i, :]) # current record
    n = length(x)
    @inbounds for j in 1:n
        if x[j] == 0
            print(pb, "\t0/0")
        elseif x[j] == 1
            print(pb, "\t1/0")
        elseif x[j] == 2
            print(pb, "\t1/1")
        else
            error("imputed genotypes can only be 0, 1, 2 but got $(x[j])")
        end
    end
    print(pb, "\n")
    nothing
end

"""
Helper function for saving a record (SNP), tracking phase information.
Here `X = X1 + X2`. 
"""
function write_snp!(pb::IOBuffer, X::Tuple, i::Int)
    X1, X2 = X[1], X[2]
    x1 = @view(X1[i, :])
    x2 = @view(X2[i, :])

    n = length(x1)
    @assert n == length(x2)
    @inbounds for j in 1:n
        if x1[j] == x2[j] == 0
            print(pb, "\t0|0")
        elseif x1[j] == 0 && x2[j] == 1
            print(pb, "\t0|1")
        elseif x1[j] == 1 && x2[j] == 0
            print(pb, "\t1|0")
        elseif x1[j] == 1 && x2[j] == 1
            print(pb, "\t1|1")
        else
            error("phased genotypes can only be 0|0, 0|1, 1|0 or 1|1 but
                got $(x1[j])|$(x2[j])")
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
    Threads.@threads for person in 1:n
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
                H1 = compressed_Hunique.CW[w1].uniqueH
                H2 = compressed_Hunique.CW[w2].uniqueH
                X[snp, person] = H1[i1, h1] + H2[i2, h2]
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
