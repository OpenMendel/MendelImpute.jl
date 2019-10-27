__precompile__()

module MendelImpute

	#needed for groupslices function
	import Base.hash
	import Base.Cartesian, Base.Cartesian.@nloops, Base.Cartesian.@nref

	using LinearAlgebra
	using StatsBase
	using ElasticArrays

	export continue_haplotype, haplopair!, haplopair, haploimpute!
	export impute!, impute2!, phase, search_breakpoint, unique_haplotypes
	export unique_haplotype_idx, groupslices, groupslices!, phase3
	export compute_redundant_haplotypes!, redundant_haplotypes
	export HaplotypeMosaicPair, HaplotypeMosaic, UniqueHaplotypeMaps
	export PeoplesRedundantHaplotypeSet, phase2, non_redundant_haplotypes

	# export UniqueHaplotypes, fast_elimination, unique_index!
	export OptimalHaplotypeSet, compute_optimal_halotype_set
	export make_refvcf_file, make_tgtvcf_file

	include("data_structures.jl")
	include("haplotyping.jl")
	include("utilities.jl")
	include("simulate_utilities.jl")

end # module
