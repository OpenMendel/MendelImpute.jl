"""
    vcf_to_jlso(vcffile, outfile, [column_major = true])

Converts haplotypes in `vcffile` into a `RefHaplotypes` object and saves it as a `jlso` format. 
It is 20~50x faster to read reference haplotypes saved in this format than `.vcf.gz` format and 
barely requires extra space to store. Internally haplotypes are stored as `BitMatrix`, so all 
haplotypes must be phased and non missing. 

# Input
* `vcffile`: file name, must end in `.vcf` or `.vcf.gz`
* `outfile`: Output file name (with or without `.jlso` extension)
* `dims`: Orientation of `H`. `2` means save haplotype vectors as columns. `1` means save as rows. 

# Output
* `hapset`: Data structure for keeping track of unique haplotypes in each window (written to `outfile`). 
"""
function vcf_to_jlso(
    vcffile::AbstractString,
    outfile::AbstractString;
    column_major::Bool = true, 
    )
    # import data
    H, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt = convert_ht(Bool, vcffile, trans=column_major, save_snp_info=true, msg = "Importing reference haplotype files...")
    
    # create `RefHaplotypes` data structure
    hapset = RefHaplotypes(H, column_major, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt)

    # save hapset to binary file using JLD2 package
    endswith(outfile, ".jlso") || (outfile = outfile * ".jlso")
    JLSO.save(outfile, :hapset => hapset, format=:julia_serialize, compression=:gzip)

    return hapset
end
