# TODO index files
function convert_gt(b::Bgen, T=Float64)
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
        X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt = 
            convert_gt(tgtfile, Float32)
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
        H, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt = convert_ht(Bool, 
            reffile, trans=true, save_snp_info=true, 
            msg="importing reference data...")
    elseif endswith(reffile, ".bgen")
        # TODO
    else
        error("Unrecognized reference file format: only VCF (ends in .vcf" * 
            " or .vcf.gz) or BGEN (ends in .bgen) files are accepted.")
    end
    return H, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt
end
