__precompile__()

module MendelImpute

    import StatsBase: sample
    
    using LinearAlgebra
    using StatsBase
    using GeneticVariation
    using VCFTools
    using SnpArrays
    using GroupSlices
    using Random
    using ProgressMeter
    using JLD2, FileIO
    using JLSO
    using Distances
    using ThreadPools
    using Lasso

    export continue_haplotype, haplopair!, haplopair, haploimpute!
    export impute!, impute_discard_phase!, search_breakpoint
    export unique_haplotype_idx
    export compute_redundant_haplotypes!, redundant_haplotypes
    export HaplotypeMosaicPair, HaplotypeMosaic, UniqueHaplotypeMaps
    export simulate_genotypes2
    export compute_optimal_halotype_pair
    export simulate_phased_genotypes
    export connect_happairs
    export phase!, phase_fast!, nearest_window_with_sufficient_typed_snps
    export haplopair_thin, haplopair_screen, haplopair_lasso

    # main functions that users are exposed to
    export phase
    export compress_haplotypes

    export OptimalHaplotypeSet, compute_optimal_halotype_set
    export make_refvcf_file, make_tgtvcf_file
    export simulate_uniform_haplotypes, simulate_markov_haplotypes, simulate_genotypes
    export unphase, compress_vcf_to_gz
    export extract_marker_info

    # imputation
    export impute_typed_only!, impute_untyped!
    export update_marker_position!

    include("compress.jl")
    include("data_structures.jl")
    include("phasing.jl")
    include("haplotype_pair.jl")
    include("haplotype_pair_screen.jl")
    include("haplotype_thinning.jl")
    include("haplotype_lasso.jl")
    include("simulate_utilities.jl")
    include("dynamic_programming.jl")
    include("impute.jl")
    include("breakpoints.jl")

end # module
