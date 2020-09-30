
# Run MendelImpute as script

If you don't want to run `MendelImpute.jl` in a Julia session (e.g. you want to run batch jobs on a cluster), you can do so by putting the code above in a Julia file. For example, in order to run with 8 threads, create a file called `impute.jl` which contains:

```julia
# place these code in a file called impute.jl
using MendelImpute, VCFTools, LinearAlgebra

# setup code goes here
reffile = ARGS[1]       # first command line argument
tgtfile = ARGS[2]       # second command line argument
outfile = ARGS[3]       # third command line argument
BLAS.set_num_threads(1) # set BLAS threads to 1 (see performance gotchas)

# run MendelImpute with default options
phase(tgtfile, reffile, outfile)
```

Then in the terminal/command-prompt, you can do
```
export JULIA_NUM_THREADS=8
julia impute.jl ref.jlso target.vcf.gz output.vcf.gz
```
