# INSTRUCTIONS: Start a Julia session on Hoffman and run code below chunk-by-chunk

target_data = ["data1.vcf.gz", "data4.vcf.gz"]
memory = 47 
missingprop = 0.1
samples = 1000

# filter/mask/unphase data. 
# THIS CHUNK HAS TO RUN BEFORE RUNNING ANY CHUNKS BELOW
cd("/u/home/b/biona001/haplotype_comparisons/data")
for data in target_data
    # if unphased target genotype file exist already, move to next step
    if isfile(data * "./tgt_masked_unphased_" * data)
        continue 
    end

    open("filter.sh", "w") do io
        println(io, "#!/bin/bash")
        println(io, "#\$ -cwd")
        println(io, "# error = Merged with joblog")
        println(io, "#\$ -o joblog.\$JOB_ID")
        println(io, "#\$ -j y")
        println(io, "#\$ -l arch=intel-X5650,exclusive,h_rt=24:00:00,h_data=$(memory)G")
        println(io, "# Email address to notify")
        println(io, "#\$ -M \$USER@mail")
        println(io, "# Notify when")
        println(io, "#\$ -m a")
        println(io)
        println(io, "# load the job environment:")
        println(io, ". /u/local/Modules/default/init/modules.sh")
        println(io, "module load julia/1.2.0")
        println(io, "module load R/3.5.1")
        println(io, "module load java/1.8.0_111")
        println(io)
        println(io, "# filter/mask data")
        println(io, "julia ./filter_and_mask.jl $data $samples")
    end

    # submit job
    run(`qsub filter.sh`)
    rm("filter.sh", force=true)
    sleep(2)
end

# AFTER above jobs finish, run prephasing using beagle 4.1. 
# THIS CHUNK HAS TO FINISH RUNNING BEFORE MINIMAC4 AND BEAGLE 5 CHUNKS CAN RUN 
for data in target_data
    open("prephase.sh", "w") do io
        println(io, "#!/bin/bash")
        println(io, "#\$ -cwd")
        println(io, "# error = Merged with joblog")
        println(io, "#\$ -o joblog.\$JOB_ID")
        println(io, "#\$ -j y")
        println(io, "#\$ -l arch=intel-X5650,exclusive,h_rt=24:00:00,h_data=$(memory)G")
        println(io, "# Email address to notify")
        println(io, "#\$ -M \$USER@mail")
        println(io, "# Notify when")
        println(io, "#\$ -m a")
        println(io)
        println(io, "# load the job environment:")
        println(io, ". /u/local/Modules/default/init/modules.sh")
        println(io, "module load julia/1.2.0")
        println(io, "module load R/3.5.1")
        println(io, "module load java/1.8.0_111")
        println(io)
        println(io, "# run prephasing using beagle 4.1")
        println(io, "java -Xss5m -Xmx$(memory)g -jar beagle4.1.jar gt=./tgt_masked_unphased_$(data) ref=./ref_$(data) niterations=0 out=./tgt_masked_phased_$(data)")
    end
    # submit job
    run(`qsub prephase.sh`)
    rm("prephase.sh", force=true)
    sleep(2)
end

# run MendelImpute (fast version and dynamic programming version)
widths = [400; 800; 1600]
for data in target_data, width in widths
    # fast version
    open("mendel$width.sh", "w") do io
        println(io, "#!/bin/bash")
        println(io, "#\$ -cwd")
        println(io, "# error = Merged with joblog")
        println(io, "#\$ -o joblog.\$JOB_ID")
        println(io, "#\$ -j y")
        println(io, "#\$ -l arch=intel-X5650,exclusive,h_rt=24:00:00,h_data=$(memory)G")
        println(io, "# Email address to notify")
        println(io, "#\$ -M \$USER@mail")
        println(io, "# Notify when")
        println(io, "#\$ -m a")
        println(io)
        println(io, "# load the job environment:")
        println(io, ". /u/local/Modules/default/init/modules.sh")
        println(io, "module load julia/1.2.0")
        println(io, "module load R/3.5.1")
        println(io, "module load java/1.8.0_111")
        println(io)
        println(io, "# run MendelImpute (fast)")
        println(io, "julia mendel_fast.jl $data $width")
    end
    # submit job
    run(`qsub mendel$width.sh`)
    rm("mendel$width.sh", force=true)
    sleep(2)

    # dynamic programming version
    open("dp_mendel$width.sh", "w") do io
        println(io, "#!/bin/bash")
        println(io, "#\$ -cwd")
        println(io, "# error = Merged with joblog")
        println(io, "#\$ -o joblog.\$JOB_ID")
        println(io, "#\$ -j y")
        println(io, "#\$ -l arch=intel-X5650,exclusive,h_rt=24:00:00,h_data=$(memory)G")
        println(io, "# Email address to notify")
        println(io, "#\$ -M \$USER@mail")
        println(io, "# Notify when")
        println(io, "#\$ -m a")
        println(io)
        println(io, "# load the job environment:")
        println(io, ". /u/local/Modules/default/init/modules.sh")
        println(io, "module load julia/1.2.0")
        println(io, "module load R/3.5.1")
        println(io, "module load java/1.8.0_111")
        println(io)
        println(io, "# run MendelImpute (dynamic programming)")
        println(io, "julia mendel_dp.jl $data $width")
    end
    # submit job
    run(`qsub dp_mendel$width.sh`)
    rm("dp_mendel$width.sh", force=true)
    sleep(2)
end

# run minimac4
for data in target_data
    open("minimac4.sh", "w") do io
        println(io, "#!/bin/bash")
        println(io, "#\$ -cwd")
        println(io, "# error = Merged with joblog")
        println(io, "#\$ -o joblog.\$JOB_ID")
        println(io, "#\$ -j y")
        println(io, "#\$ -l arch=intel-X5650,exclusive,h_rt=24:00:00,h_data=$(memory)G")
        println(io, "# Email address to notify")
        println(io, "#\$ -M \$USER@mail")
        println(io, "# Notify when")
        println(io, "#\$ -m a")
        println(io)
        println(io, "# load the job environment:")
        println(io, ". /u/local/Modules/default/init/modules.sh")
        println(io, "module load julia/1.2.0")
        println(io, "module load R/3.5.1")
        println(io, "module load java/1.8.0_111")
        println(io)
        println(io, "# Impute using minimac4")
        println(io, "java -Xmx$(memory)g -jar beagle5.0.jar gt=./tgt_masked_unphased_$(data) ref=./ref_$(data) out=./beagle_imputed_$(data)")
    end
    # submit job
    run(`qsub minimac4.sh`)
    rm("minimac4.sh", force=true)
    sleep(2)
end

# run beagle 5.0 
for data in target_data
    open("beagle.sh", "w") do io
        println(io, "#!/bin/bash")
        println(io, "#\$ -cwd")
        println(io, "# error = Merged with joblog")
        println(io, "#\$ -o joblog.\$JOB_ID")
        println(io, "#\$ -j y")
        println(io, "#\$ -l arch=intel-X5650,exclusive,h_rt=24:00:00,h_data=$(memory)G")
        println(io, "# Email address to notify")
        println(io, "#\$ -M \$USER@mail")
        println(io, "# Notify when")
        println(io, "#\$ -m a")
        println(io)
        println(io, "# load the job environment:")
        println(io, ". /u/local/Modules/default/init/modules.sh")
        println(io, "module load julia/1.2.0")
        println(io, "module load R/3.5.1")
        println(io, "module load java/1.8.0_111")
        println(io)
        println(io, "# Impute using beagle 5.0")
        println(io, "java -Xmx$(memory)g -jar beagle5.0.jar gt=./tgt_masked_unphased_$(data) ref=./ref_$(data) out=./beagle_imputed_$(data)")
    end
    # submit job
    run(`qsub beagle.sh`)
    rm("beagle.sh", force=true)
    sleep(2)
end



