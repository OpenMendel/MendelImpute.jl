"""
Runs Impute 5 on chromosome 20 over 4 chunks.
Each chunk imputes 20 million SNPs using 2 or 3 threads.
Each script runs a different chunk and are initialized simultaneously.
"""
function impute5_HRC_chr20(
    reffile::AbstractString,
    tgtfile::AbstractString,
    )
    chr = 20                # chromosome 20
    allale_start = 1        # first SNP is 60309
    allele_end = 70000000   # last SNP is 69960292
    chunk_size = 20000000   # each chunk imputes 20cM
    nchunks = 4             # chunk 1 imputes first 20cM region, chunk 2 imputes next 20cM ...etc
    cthreads = [3, 3, 2, 2] # chunk 1 uses 3 threads, chunk 2 uses 3 threads...etc

    # start 4 processes, each imputing different regions, using 2-3 threads
    Threads.@threads for i in 1:nchunks
        cur_start = (i - 1) * chunk_size + 1
        cur_end = i == nchunks ? allele_end : chunk_size * i
        cur_threads = cthreads[i]
        run(`nohup /usr/bin/time /home/biona001/impute5_v1.1.4/impute5_1.1.4_static
            --h $reffile
            --g $tgtfile
            --r $chr:$cur_start-$cur_end
            --o impute5.chr$chr.result.chunk$i.vcf.gz
            --threads $cur_threads`)
    end
end

"""
Runs Impute 5 on chromosome 10 over 14 chunks.
Each chunk imputes 10 million SNPs using 1 thread.
Up to 8 chunks can be ran simultaneously (10 concurrent threads throw memory error).

Note using 20cM chunks (as recommended) would result in 7 chunks, 
"""
function impute5_HRC_chr10(
    reffile::AbstractString,
    tgtfile::AbstractString,
    )
    chr = 10                # chromosome 10
    allale_start = 1        # first SNP is 94316
    allele_end = 140000000  # last SNP is 135523061
    chunk_size = 10000000   # each chunk imputes 10cM (running 20cM chunks with 7 chunks throws Our Of Memory Error)
    nchunks = 14            # chunk 1 imputes first 10cM region, chunk 2 imputes next 10cM ...etc
    cthreads = ones(Int, nchunks)  # each chunk runs 1 threads

    # run each chunk with 1 thread
    Threads.@threads for i in 1:nchunks
        cur_start = (i - 1) * chunk_size + 1
        cur_end = i == nchunks ? allele_end : chunk_size * i
        cur_threads = cthreads[i]
        run(`/home/biona001/impute5_v1.1.4/impute5_1.1.4_static
            --h $reffile
            --g $tgtfile
            --r $chr:$cur_start-$cur_end
            --o impute5.chr$chr.result.chunk$i.vcf.gz
            --threads $cur_threads`)
    end
end

chr = parse(Int, ARGS[1])
reffile = ARGS[2]
tgtfile = ARGS[3]

if Threads.nthreads() == 1
    error("start Julia with >1 thread for this script to work!")
end

if chr == 20
    _, t, _, _, _ = @timed impute5_HRC_chr20(reffile, tgtfile)
elseif chr == 10
    _, t, _, _, _ = @timed impute5_HRC_chr10(reffile, tgtfile)
end

println("Total elapsed time = $t seconds")

