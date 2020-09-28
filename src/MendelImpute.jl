__precompile__()

module MendelImpute

    using LinearAlgebra
    using GeneticVariation
    using VCFTools
    using SnpArrays
    using GroupSlices
    using Random
    using ProgressMeter
    using JLD2, FileIO, JLSO
    using Distances
    using LazyArrays

    # main functions that users are exposed to
    export phase               # main function for imputation and phasing
    export compress_haplotypes # for compressing haplotype panels
    export paint, composition, unique_populations # chromosome painting and admix proportions
    export convert_compressed  # for ultra-compression outputs

    include("compress.jl")
    include("data_structures.jl")
    include("phasing.jl")
    include("haplotype_pair.jl")
    include("haplotype_rescreen.jl")
    include("haplotype_thinning.jl")
    include("haplotype_stepscreen.jl")
    include("dynamic_programming.jl")
    include("impute.jl")
    include("breakpoints.jl")
    include("intersect.jl")
    include("painting.jl")
    include("ultra_compress.jl")

    # test data directory
    datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)    
end # module
