module MendelImpute

    using LinearAlgebra
    using VCFTools
    using SnpArrays
    using BGEN
    using GroupSlices
    using Random
    using ProgressMeter
    using JLD2, FileIO, JLSO
    using Distances
    using LazyArrays
    using VariantCallFormat
    using CSV, DataFrames

    # functions for phasing/imputation/data-compression
    export phase               # main function for imputation and phasing
    export compress_haplotypes # for compressing haplotype panels
    export convert_compressed  # for ultra-compression outputs

    # functions for ancestry inference
    export admixture_global
    export paint, composition, unique_populations 
    export thousand_genome_population_to_superpopulation
    export thousand_genome_samples_to_population
    export thousand_genome_samples_to_super_population

    include("compress.jl")
    include("data_structures.jl")
    include("phase.jl")
    include("haplotype_pair.jl")
    include("haplotype_rescreen.jl")
    include("haplotype_thinning.jl")
    include("haplotype_stepscreen.jl")
    include("dynamic_programming.jl")
    include("impute.jl")
    include("breakpoints.jl")
    include("intersect.jl")
    include("ancestry.jl")
    include("ultra_compress.jl")
    include("read.jl")

    # test data directory
    datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)    
end # module
