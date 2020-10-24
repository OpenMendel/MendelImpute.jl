using ADMIXTURE

println("running chr18.uniqueSNPs.bed with 26 populations and 10 threads")

P, Q = admixture("chr18.uniqueSNPs.bed", 26, j=10)
