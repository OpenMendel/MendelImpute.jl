using VCFTools
using MendelImpute
using Random

"""
    filter_and_mask(data::String, samples::Int)

Creates reference haplotypes and (unphased) target genotype files from `data`. 

# Inputs
`data`: The full (phased) data simulated by msprime.
`samples`: Number of samples (genotypes) desired in target file. Remaining haplotypes will become the reference panel
"""
function filter_and_mask(data::String, samples::String)
    missingprop = 0.1
    n = nsamples(data)
    p = nrecords(data)
    samples = convert
    samples > n && error("requested samples exceed total number of genotypes in $data.")

    # output filenames (tgt_data1.vcf.gz, ref_data1.vcf.gz, and tgt_masked_data1.vcf.gz)
    tgt = "./tgt_" * data
    ref = "./ref_" * data
    tgt_mask = "./tgt_masked_" * data
    tgt_mask_unphase = "./tgt_masked_unphased_" * data

    # compute target and reference index
    tgt_index = falses(n)
    tgt_index[1:samples] .= true
    ref_index = .!tgt_index
    record_index = 1:p # save all records (SNPs) 

    # generate masking matrix with `missingprop`% of trues (true = convert to missing)
    Random.seed!(2020)
    masks = falses(p, samples)
    for j in 1:samples, i in 1:p
        rand() < missingprop && (masks[i, j] = true)
    end

    # create outputs 
    VCFTools.filter(data, record_index, tgt_index, des = tgt)
    VCFTools.filter(data, record_index, ref_index, des = ref)
    mask_gt(tgt, masks, des=tgt_mask)

    # finally, unphase the target data
    unphase(tgt_mask, outfile=tgt_mask_unphase)
end

data = ARGS[1]
samples = parse(Int, ARGS[2])
filter_and_mask(data, samples)
