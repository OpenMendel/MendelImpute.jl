###### This file is part of the MendelImpute.jl package.
###### These are wrapper functions for importing VCF/PLINK/BGEN target
###### and reference files. 

function convert_gt(b::Bgen, T=Float32)
    n = n_samples(b)
    p = n_variants(b)

    # return arrays
    G = Matrix{T}(undef, p, n)
    Gchr = Vector{String}(undef, p)
    Gpos = Vector{String}(undef, p)
    GsnpID = Vector{String}(undef, p)
    Gref = Vector{String}(undef, p)
    Galt = Vector{String}(undef, p)

    # loop over each variant
    i = 1
    for v in iterator(b; from_bgen_start=true)
        dose = minor_allele_dosage!(b, v; T=T)
        copyto!(@view(G[i, :]), dose)
        Gchr[i], pos[i], GsnpID[i], Gref[i], Galt[i] =
            chrom(v), pos(v), rsid(v), major_allele(v), minor_allele(v)
        i += 1
        clear!(v)
    end
    return G, chr, pos, snpID, ref, alt
end

function convert_ht(b::Bgen)
    n = 2n_samples(b)
    p = n_variants(b)

    # return arrays
    H = BitMatrix(undef, p, n)
    Hchr = Vector{String}(undef, p)
    Hpos = Vector{String}(undef, p)
    HsnpID = Vector{String}(undef, p)
    Href = Vector{String}(undef, p)
    Halt = Vector{String}(undef, p)

    # loop over each variant
    i = 1
    for v in iterator(b; from_bgen_start=true)
        dose = probabilities!(b, v; T=Bool)
        # copyto!(@view(G[i, :]), dose)
        # Hchr[i], Hpos[i], HsnpID[i], Href[i], Halt[i] =
        #     chrom(v), pos(v), rsid(v), major_allele(v), minor_allele(v)
        i += 1
        clear!(v)
    end
    return H, chr, pos, snpID, ref, alt
end

isplink(tgtfile::AbstractString) = isfile(tgtfile * ".bed") && 
                                   isfile(tgtfile * ".fam") && 
                                   isfile(tgtfile * ".bim")

function import_target(tgtfile::AbstractString, dosage=false)
    if (endswith(tgtfile, ".vcf") || endswith(tgtfile, ".vcf.gz")) && !dosage
        X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt = 
            VCFTools.convert_gt(UInt8, tgtfile, trans=true, 
            save_snp_info=true, msg = "Importing genotype file...")
    elseif (endswith(tgtfile, ".vcf") || endswith(tgtfile, ".vcf.gz")) && dosage
        X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt = 
            VCFTools.convert_ds(Float32, tgtfile, trans=true, 
            save_snp_info=true, msg = "Importing genotype file as dosages...")
    elseif endswith(tgtfile, ".bgen")
        samplefile = tgtfile[1:end-5] * ".sample"
        isfile(samplefile) || error("sample file $samplefile not found!")
        indexfile = isfile(tgtfile * ".bgi") ? tgtfile * ".bgi" : nothing
        bgen = Bgen(tgtfile; sample_path=samplefile, idx_path=indexfile)
        X, X_chr, X_pos, X_ids, X_ref, X_alt = convert_gt(bgen, Float32)
        X_sampleID = BGEN.get_samples(samplefile, n_samples(b))
    elseif isplink(tgtfile)
        dosage && error("PLINK files detected but dosage = true!")
        # convert SnpArray data to matrix.
        X_snpdata = SnpArrays.SnpData(tgtfile)
        X = convert(Matrix{Union{UInt8, Missing}}, Transpose(X_snpdata.snparray))
        X[findall(isone, X)] .= missing     # 0x01 encodes missing
        X[findall(x -> x === 0x02, X)] .= 1 # 0x02 is 1
        X[findall(x -> x === 0x03, X)] .= 2 # 0x03 is 2
        # get other relevant information
        X_sampleID = X_snpdata.person_info[!, :iid]
        X_chr = X_snpdata.snp_info[!, :chromosome]
        X_pos = X_snpdata.snp_info[!, :position]
        X_ids = X_snpdata.snp_info[!, :snpid]
        X_ref = X_snpdata.snp_info[!, :allele1]
        X_alt = X_snpdata.snp_info[!, :allele2]
    else
        error("Unrecognized target file format: target file can only be VCF" *
            " files (ends in .vcf or .vcf.gz), BGEN (ends in .bgen) or PLINK" *
            " (do not include.bim/bed/fam) and all trio must exist in 1 directory)")
    end
    return X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt
end

function import_compressed_reference(reffile::AbstractString)
    println("Importing reference haplotype data...")
    flush(stdout)
    
    if endswith(reffile, ".jlso")
        compressed_Hunique = read_jlso(reffile)
    elseif endswith(reffile, ".vcf") || endswith(reffile, ".vcf.gz")
        error("reference panel is VCF format: please first compress reference" *
            " file to .jlso format using the compress_haplotypes() function.")
    elseif endswith(reffile, ".bgen")
        error("reference panel is BGEN format: please first compress reference" *
            " file to .jlso format using the compress_haplotypes() function.")
    elseif endswith(reffile, ".jld2")
        @load reffile compressed_Hunique 
    else
        error("Unrecognized reference file format: only VCF (ends in .vcf" * 
            " or .vcf.gz), BGEN (ends in .bgen) `.jlso` files are accepted.")
    end
    return compressed_Hunique
end

function import_reference(reffile::AbstractString)
    if endswith(reffile, ".vcf") || endswith(reffile, ".vcf.gz")
        H, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt = 
            VCFTools.convert_ht(Bool, reffile, trans=true, save_snp_info=true, 
            msg="importing reference data...")
    elseif endswith(reffile, ".bgen")
        samplefile = reffile[1:end-5] * ".sample"
        isfile(samplefile) || error("sample file $samplefile not found!")
        indexfile = isfile(reffile * ".bgi") ? reffile * ".bgi" : nothing
        bgen = Bgen(reffile; sample_path=samplefile, idx_path=indexfile)
        H, H_chr, H_pos, H_ids, H_ref, H_alt = convert_ht(bgen)
        H_sampleID = BGEN.get_samples(samplefile, n_samples(b))
    else
        error("Unrecognized reference file format: only VCF (ends in .vcf" * 
            " or .vcf.gz) or BGEN (ends in .bgen) files are accepted.")
    end
    return H, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt
end
