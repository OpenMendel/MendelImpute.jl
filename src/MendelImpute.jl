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
    using JLSO
    using Distances
    using ThreadPools

    export continue_haplotype, haplopair!, haplopair, haploimpute!
    export impute!, impute_discard_phase!, search_breakpoint
    export unique_haplotype_idx
    export compute_redundant_haplotypes!, redundant_haplotypes
    export HaplotypeMosaicPair, HaplotypeMosaic, UniqueHaplotypeMaps
    export simulate_genotypes2
    export compute_optimal_halotype_pair
    export simulate_phased_genotypes
    export connect_happairs
    export phase!, phase_fast!, phase_sample!
    export haplopair_screen!, haplopair_stepscreen!

    # main functions that users are exposed to
    export phase
    export compress_haplotypes
    export nhaplotypes, windows, count_haplotypes_per_window
    export avg_haplotypes_per_window, nchunks, max_haplotypes_per_window

    export compute_optimal_haplotypes!
    export make_refvcf_file, make_tgtvcf_file
    export simulate_uniform_haplotypes, simulate_markov_haplotypes
    export simulate_genotypes
    export unphase, compress_vcf_to_gz
    export extract_marker_info

    # imputation
    export impute_typed_only!, impute_untyped!
    export update_marker_position!

    include("compress.jl")
    include("data_structures.jl")
    include("phasing.jl")
    include("haplotype_pair.jl")
    include("haplotype_rescreen.jl")
    include("haplotype_thinning.jl")
    include("haplotype_stepscreen.jl")
    include("simulate_utilities.jl")
    include("dynamic_programming.jl")
    include("impute.jl")
    include("breakpoints.jl")
    include("intersect.jl")

end # module
