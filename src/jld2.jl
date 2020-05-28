"""
    save_jld2(vcffile, outfile, [width], [dims], [flankwidth])

Converts haplotypes in `vcffile` into a `RefHaplotypes` object and saves it as a `jld2` format. 
It is ~500x faster to read reference haplotypes saved in this format than `.vcf.gz` format. 
All haplotypes must be phased and non missing. 

# Input
* `vcffile`: file name, must end in `.vcf` or `.vcf.gz`
* `outfile`: Output file name (with or without `.jld2` extension)
* `dims`: Orientation of `H`. `2` means save haplotype vectors as columns. `1` means save as rows. 

# Output
* `hapset`: Data structure for keeping track of unique haplotypes in each window (written to `outfile`). 
"""
function save_jld2(
    vcffile::AbstractString,
    outfile::AbstractString;
    column_major::Bool = true, 
    )
    # import data
    H, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt = convert_ht(Bool, vcffile, trans=column_major, save_snp_info=true, msg = "Importing reference haplotype files...")
    
    # create `RefHaplotypes` data structure
    hapset = RefHaplotypes(H, column_major, H_sampleID, H_chr, H_pos, H_ids, H_ref, H_alt)

    # save hapset to binary file using JLD2 package
    endswith(outfile, ".jld2") || (outfile = outfile * ".jld2")
    @save outfile hapset

    return hapset
end
