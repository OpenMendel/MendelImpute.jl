###### This file is part of the MendelImpute.jl package.
###### These are wrapper functions for importing VCF/PLINK/BGEN target
###### and reference files. 

"""
    convert_gt(b::Bgen, T=Float32)

Imports dosage information and chr/pos/snpID/ref/alt into numeric arrays.

# Input
- `b`: a `Bgen` object
- `T`: Type for genotype array

# Output
- `G`: a `p × n` matrix of type `T`. Each column is a genotype
- `Gchr`: Vector of `String`s holding chromosome number for each variant
- `Gpos`: Vector of `Int` holding each variant's position
- `GsnpID`: Vector of `String`s holding variant ID for each variant
- `Gref`: Vector of `String`s holding reference allele for each variant
- `Galt`: Vector of `String`s holding alterante allele for each variant
"""
function convert_gt(b::Bgen, T=Float32)
    n = n_samples(b)
    p = n_variants(b)

    # return arrays
    G = Matrix{T}(undef, p, n)
    Gchr = Vector{String}(undef, p)
    Gpos = Vector{Int}(undef, p)
    GsnpID = Vector{String}(undef, p)
    Gref = Vector{String}(undef, p)
    Galt = Vector{String}(undef, p)

    # loop over each variant
    i = 1
    for v in iterator(b; from_bgen_start=true)
        dose = minor_allele_dosage!(b, v; T=T)
        copyto!(@view(G[i, :]), dose)
        Gchr[i], Gpos[i], GsnpID[i], Gref[i], Galt[i] =
            chrom(v), pos(v), rsid(v), major_allele(v), minor_allele(v)
        i += 1
        clear!(v)
    end
    return G, Gchr, Gpos, GsnpID, Gref, Galt
end

"""
    convert_ht(b::Bgen)

Import phased haplotypes as a `BitMatrix`, and store chr/pos/snpID/ref/alt.

# Input
- `b`: a `Bgen` object. Each variant must be phased and samples must be diploid

# Output
- `H`: a `p × 2n` matrix of type `T`. Each column is a haplotype. 
- `Hchr`: Vector of `String`s holding chromosome number for each variant
- `Hpos`: Vector of `Int` holding each variant's position
- `HsnpID`: Vector of `String`s holding variant ID for each variant
- `Href`: Vector of `String`s holding reference allele for each variant
- `Halt`: Vector of `String`s holding alterante allele for each variant
"""
function convert_ht(b::Bgen)
    n = 2n_samples(b)
    p = n_variants(b)

    # return arrays
    H = BitMatrix(undef, p, n)
    Hchr = Vector{String}(undef, p)
    Hpos = Vector{Int}(undef, p)
    HsnpID = Vector{String}(undef, p)
    Href = Vector{String}(undef, p)
    Halt = Vector{String}(undef, p)

    # loop over each variant
    i = 1
    for v in iterator(b; from_bgen_start=true)
        dose = probabilities!(b, v)
        phased(v) || error("variant $(rsid(v)) at position $(pos(v)) not phased!")
        for j in 1:n_samples(b)
            Hi = @view(dose[:, j])
            H[i, 2j - 1] = read_haplotype1(Hi)
            H[i, 2j] = read_haplotype2(Hi)
        end
        Hchr[i], Hpos[i], HsnpID[i], Href[i], Halt[i] =
            chrom(v), pos(v), rsid(v), major_allele(v), minor_allele(v)
        i += 1
        clear!(v)
    end
    return H, Hchr, Hpos, HsnpID, Href, Halt
end

read_haplotype1(Hi::AbstractVector) = Hi[2] ≥ 0.5 ? true : false
read_haplotype2(Hi::AbstractVector) = Hi[4] ≥ 0.5 ? true : false

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
        # isfile(samplefile) || error("sample file $samplefile not found!")
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

"""
    read_jlso(file::AbstractString)

Imports a `.jlso`-compressed reference haplotype panel.
"""
function read_jlso(reffile::AbstractString)
    loaded = JLSO.load(reffile)
    return loaded[:compressed_Hunique]
end
