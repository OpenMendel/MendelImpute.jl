__precompile__()

module MendelImpute

    import StatsBase: sample
    
    using LinearAlgebra
    using StatsBase
    using ElasticArrays
    using GeneticVariation
    using VCFTools
    using GroupSlices

    export continue_haplotype, haplopair!, haplopair, haploimpute!
    export impute!, impute2!, phase, search_breakpoint, unique_haplotypes
    export unique_haplotype_idx
    export compute_redundant_haplotypes!, redundant_haplotypes
    export HaplotypeMosaicPair, HaplotypeMosaic, UniqueHaplotypeMaps
    export PeoplesRedundantHaplotypeSet, phase2, non_redundant_haplotypes
    export simulate_genotypes2
    export set_flip!
    export phase_prephased
    export compute_optimal_halotype_set_prephased
    export simulate_phased_genotypes

    # export UniqueHaplotypes, fast_elimination, unique_index!
    export OptimalHaplotypeSet, compute_optimal_halotype_set
    export make_refvcf_file, make_tgtvcf_file
    export simulate_uniform_haplotypes, simulate_markov_haplotypes, simulate_genotypes

    include("data_structures.jl")
    include("haplotyping.jl")
    include("utilities.jl")
    include("simulate_utilities.jl")
    include("dynamic_programming.jl")

end # module
