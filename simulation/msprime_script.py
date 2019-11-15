import msprime, sys

# parameters
sample_size = int(sys.argv[1])
effective_population_size = int(sys.argv[2])
sequence_length = int(sys.argv[3])
recombination_rate=float(sys.argv[4])
mutation_rate=float(sys.argv[5])
seed = int(sys.argv[6])

# run the simulation
tree_sequence = msprime.simulate(sample_size=sample_size, Ne=effective_population_size, length=sequence_length, recombination_rate=recombination_rate, mutation_rate=mutation_rate, random_seed=seed)

# print results
with sys.stdout as vcffile:
	tree_sequence.write_vcf(vcffile, 2) # 2 is for diploid
