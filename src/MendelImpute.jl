__precompile__()

module MendelImpute

	#needed for groupslices function
	import Base.hash
	import Base.Cartesian, Base.Cartesian.@nloops, Base.Cartesian.@nref

	using LinearAlgebra
	using StatsBase
	using ElasticArrays

	export continue_haplotype, haplopair!, haplopair, haploimpute!
	export impute!, phase, search_breakpoint, unique_haplotypes
	export unique_haplotype_idx, groupslices, groupslices!
	export compute_redundant_haplotypes!, redundant_haplotypes
	export HaplotypeMosaicPair, HaplotypeMosaic, UniqueHaplotypeMaps, PeoplesRedundantHaplotypeSet

	include("data_structures.jl")
	include("haplotyping.jl")

end # module
