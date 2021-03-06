using VCFTools
using MendelImpute
using GeneticVariation

function run(data::String, width::Int)
    tgtfile = "./tgt_masked_unphased_" * data
    reffile = "./ref_" * data
    outfile = "./mendel_imputed_" * data
    phase(tgtfile, reffile, outfile = outfile, width = width, fast_method=true)
end

data  = ARGS[1]
width = ARGS[2]
run(data, width)
