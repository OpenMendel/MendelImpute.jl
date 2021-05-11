
# Performance gotchas

## Gotcha 1: Run MendelImpute in parallel

To run `MendelImpute.jl` in parallel,
1. To use 4 threads, execute `export JULIA_NUM_THREADS=4` **before** starting Julia. 
2. Verify the Julia session is running in parallel by executing `Threads.nthreads()` in Julia
3. Set the number of BLAS threads to be 1 by `using LinearAlgebra; BLAS.set_num_threads(1)`. This avoids [oversubscription](https://ieeexplore.ieee.org/document/5470434). 

!!! note

    We recommend number of threads equal to the number of physical CPU cores on your machine. **Number of Julia threads should never exceed number of physical CPU cores!!** Hyperthreading is valuable for I/O operations (in our experience), but not for linear algebra routines used throughout MendelImpute. 

## Gotcha 2: Compressing haplotype panels is slow

Currently it is recommended to build a new compressed reference haplotype panel for every new set of typed SNPs (although this is not strictly required). The compression routine is slow because reading raw VCF files is slow. Thus, it is highly advised that one try to use the same set of typed SNPs as much as possible. 

We are actively developing a new set of functions in [SnpArrays.jl](https://github.com/OpenMendel/SnpArrays.jl) to alleviate this problem. Since SnpArrays use memory mapping, read times can be improved dramatically. 

## Gotcha 3: `max_d` too high (or too low)

When you compress the haplotype panels into a `.jlso` format, you specified `max_d` which is the maximum number of unique haplotypes per window. We generally recommend using `max_d = 1000`, BUT 1000 may be too small if you use a reference panel larger than HRC. In that case, you can try larger `max_d`, which will improve error rate. 

### Symptoms for `max_d` too high:

`Computing optimal haplotypes` is too slow. In particular, the timing for `haplopair search` is too high. 

### Symptoms for `max_d` too low:

Too few typed SNPs per window indicates `max_d` is set too low. You can calculate the number of typed SNPs per window by dividing the total number of SNPs in the target file by the total windows (a number that will be output after every run). Ideally you want an average of 400 typed SNPs per window, but something as low as 50 still works. Something like 10~20 is too low. 

### I really want to use a high `max_d`

A high `max_d` generally improve error, so it is understandable you want to do so. If a high `max_d` value runs too slow, try setting `stepwise = 100` and `max_haplotypes` to a number that is close to 1000. This avoids searching the global minimizer of the least-squares problem for windows that have more than `max_haplotypes` number of unique haplotypes. Setting `thinning_factor` instead of `stepwise` have a similar effect. Details for these 2 heuristic searches are explained in the appendix of our paper. 

## Gotcha 4: Do you have enough memory (RAM)?

While MendelImpute uses the least RAM compared to competing softwares (as of 2020), it is still possible for large imputation problems to consume all available RAM. If this happens, Julia will first try to use `swap` before crashing (until all of `swap` is consumed). Monitor your RAM usage constantly to make sure this doesn't happen. On Mac/Linux machines, the `top` or `htop` command will monitor this information. Alternatively, the `/usr/bin/time` command will automatically records max RAM usage for job and whether any `swap` had been performed. 

### Rough estimate for amount of RAM needed

There are 4 things that require lots of memory:
+ The target genotype matrix $\mathbf{X}_{n \times p}$ requires $n \times p \times 8$ bits. If $\mathbf{X}$ is dosage data, then you need instead $n \times p \times 32$ bits
+ The matrix $\mathbf{M}_{d \times d}$ requires $c \times d \times d \times 32$ bits, where $c$ is the number of parallel threads used and $d$ is the number specified in the [compress_haplotypes](https://openmendel.github.io/MendelImpute.jl/dev/man/api/#MendelImpute.compress_haplotypes) function.
+ The matrix $\mathbf{N}_{n \times d}$ requires $c \times n \times d \times 32$ bits, where $c$ is the number of parallel threads used and $d$ is the number specified in the [compress_haplotypes](https://openmendel.github.io/MendelImpute.jl/dev/man/api/#MendelImpute.compress_haplotypes) function.
+ The compressed reference haplotype panel produced by the [compress_haplotypes](https://openmendel.github.io/MendelImpute.jl/dev/man/api/#MendelImpute.compress_haplotypes) function. This typically requires about $3r$ gigabytes of RAM where $r$ is your panel's size in `.vcf.gz`. 

If you do not have the above issues and your code is still running slow, file an issue on GitHub and we will take a look at it ASAP. 
