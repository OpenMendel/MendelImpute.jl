"""
    phase(tgtfile, reffile; [outfile], [width], [flankwidth], [fast_method])

Main function of MendelImpute program. Phasing (haplotying) of `tgtfile` from a pool 
of haplotypes `reffile` by sliding windows and saves result in `outfile`. 

# Input
- `tgtfile`: VCF or PLINK files. VCF files should end in `.vcf` or `.vcf.gz`. PLINK files should exclude `.bim/.bed/.fam` suffixes but the trio must all be present in the directory. 
- `reffile`: VCF or compressed Julia binary files. VCF files should end in `.vcf` or `.vcf.gz`. Acceptable Julia binary formats includes `.jld2` (fastest read time) and `.jlso` (smallest file size).

# Optional Inputs
- `outfile`: output filename ending in `.vcf.gz` or `.vcf`. Output genotypes will have no missing data.
- `impute`: If `true`, untyped SNPs will be imputed, otherwise only missing snps in `tgtfile` will be imputed.  (default `false`)
- `width`: number of SNPs (markers) in each haplotype window. (default `2048`)
- `flankwidth`: Number of SNPs flanking the sliding window (defaults to 10% of `width`)
- `min_typed_snps`: Number of typed SNPs required in each window. Below this threshold means optimal haplotype pair will be copied from the nearest window with enough typed SNPs
"""
function phase(
    tgtfile::AbstractString,
    reffile::AbstractString;
    outfile::AbstractString = "imputed." * tgtfile,
    impute::Bool = true,
    width::Int = 2048,
    min_typed_snps = 50, 
    )

    # decide how to partition the data based on available memory 
    # people = nsamples(tgtfile)
    # haplotypes = 2nsamples(reffile_aligned)
    # snps_per_chunk = chunk_size(people, haplotypes)
    # chunks = ceil(Int, tgt_snps / snps_per_chunk)
    # remaining_snps = tgt_snps - ((chunks - 1) * snps_per_chunk)
    # println("Running chunk $chunks / $chunks")

    # import reference data
    import_data_start = time()
    if endswith(reffile, ".jld2")
        @load reffile compressed_Hunique 
        width == compressed_Hunique.width || error("Specified width = $width does not equal $(compressed_Hunique.width) = width in .jdl2 file")
    elseif endswith(reffile, ".jlso")
        loaded = JLSO.load(reffile)
        compressed_Hunique = loaded[:compressed_Hunique]
        width == compressed_Hunique.width || error("Specified width = $width does not equal $(compressed_Hunique.width) = width in .jlso file")
    elseif endswith(reffile, ".vcf") || endswith(reffile, ".vcf.gz")
        # for VCF files, compress into jlso files and filter for unique haplotypes in each window
        @info "VCF files detected: compressing reference file to .jlso format..."
        compressed_Hunique = compress_haplotypes(reffile, tgtfile, "compressed." * reffile, width, dims=2, flankwidth = 0)
    else
        error("Unrecognized reference file format: only VCF (ends in .vcf or .vcf.gz), `.jlso`, or `.jld2` files are acceptable.")
    end

    # import genotype data
    if endswith(tgtfile, ".vcf") || endswith(tgtfile, ".vcf.gz")
        X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt = VCFTools.convert_gt(UInt8, tgtfile, trans=true, save_snp_info=true, msg = "Importing genotype file...")
    elseif isfile(tgtfile * ".bed") && isfile(tgtfile * ".fam") && isfile(tgtfile * ".bim")
        # PLINK files
        X_snpdata = SnpArrays.SnpData(tgtfile)
        X = convert(Matrix{UInt8}, X_snpdata.snparray') # transpose genotypes from rows to column 
        X_sampleID = X_snpdata.person_info[!, :iid] 
        X_chr = X_snpdata.snp_info[!, :chromosome]
        X_pos = X_snpdata.snp_info[!, :position]
        X_ids = X_snpdata.snp_info[!, :snpid]
        X_ref = X_snpdata.snp_info[!, :allele1]
        X_alt = X_snpdata.snp_info[!, :allele2]
    else
        error("Unrecognized target file format: target file can only be VCF files (ends in .vcf or .vcf.gz) or PLINK files (do not include .bim/bed/fam and all three files must exist in 1 directory)")
    end
    import_data_time = time() - import_data_start

    # declare some constants
    people = size(X, 2)
    tgt_snps = size(X, 1)
    ref_snps = length(compressed_Hunique.pos)
    windows = floor(Int, tgt_snps / width)
    quad_timers = [zeros(4) for _ in 1:Threads.nthreads()]

    #
    # compute redundant haplotype sets
    #
    calculate_happairs_start = time()
    pmeter = Progress(windows, 5, "Computing optimal haplotype pairs...")
    redundant_haplotypes = [[Tuple{Int, Int}[] for i in 1:windows] for j in 1:people]
    [[sizehint!(redundant_haplotypes[j][i], 1000) for i in 1:windows] for j in 1:people] # don't save >1000 redundant happairs
    num_unique_haps = zeros(Threads.nthreads())
    Threads.@threads for w in 1:windows
        Hw_aligned = compressed_Hunique.CW_typed[w].uniqueH
        Xw_idx_start = (w - 1) * width + 1
        Xw_idx_end = (w == windows ? length(X_pos) : w * width)
        Xw_aligned = X[Xw_idx_start:Xw_idx_end, :]

        # computational routine
        happairs, hapscore, t1, t2, t3 = haplopair(Xw_aligned, Hw_aligned)
        # happairs, hapscore, t1, t2, t3 = (size(Hw_aligned, 2) < 1000 ? haplopair(Xw_aligned, Hw_aligned) :  
        #     haplopair_thin(Xw_aligned, Hw_aligned, keep=1000))
        
        # convert happairs (which index off unique haplotypes) to indices of full haplotype pool, and find all matching happairs
        t4 = @elapsed compute_redundant_haplotypes!(redundant_haplotypes, compressed_Hunique, happairs, w)

        # record timings
        id = Threads.threadid()
        quad_timers[id][1] += t1
        quad_timers[id][2] += t2
        quad_timers[id][3] += t3
        quad_timers[id][4] += t4
        num_unique_haps[id] += size(Hw_aligned, 2)

        # update progress
        next!(pmeter)
    end
    avg_num_unique_haps = sum(num_unique_haps) / windows
    calculate_happairs_time = time() - calculate_happairs_start

    #
    # phasing (haplotyping) step
    #
    # offset = (chunks - 1) * snps_per_chunk
    phase_start = time()
    ph = [HaplotypeMosaicPair(ref_snps) for i in 1:people]
    phase!(ph, X, compressed_Hunique, redundant_haplotypes, X_pos) # phase by dynamic programming + breakpoint search
    phase_time = time() - phase_start

    #
    # impute step
    #
    impute_start = time()
    H_pos = compressed_Hunique.pos
    XtoH_idx = indexin(X_pos, H_pos) # X_pos[i] == H_pos[XtoH_idx[i]]
    if impute
        # initialize whole genotype matrix and copy known entries into it
        X_full = Matrix{Union{Missing, UInt8}}(missing, ref_snps, people)
        copyto!(@view(X_full[XtoH_idx, :]), X)

        # convert phase's starting position from X's index to H's index
        update_marker_position!(ph, XtoH_idx)

        impute!(X_full, compressed_Hunique, ph, outfile, X_sampleID, XtoH_idx=nothing) # imputes X_full and writes to file
    else
        impute!(X, compressed_Hunique, ph, outfile, X_sampleID, XtoH_idx=XtoH_idx) # imputes X (only containing typed snps) and writes to file
    end
    impute_time = time() - impute_start

    println("Total windows = $windows, averaging ~ $(round(Int, avg_num_unique_haps)) unique haplotypes per window.\n")
    println("Timings: ")
    println("    Data import                     = ", round(import_data_time, sigdigits=6), " seconds")
    println("    Computing haplotype pair        = ", round(calculate_happairs_time, sigdigits=6), " seconds")
    println("        BLAS3 mul! to get M and N      = ", round(quad_timers[1][1], sigdigits=6), " seconds (on thread 1)")
    println("        haplopair search               = ", round(quad_timers[1][2], sigdigits=6), " seconds (on thread 1)")
    println("        supplying constant terms       = ", round(quad_timers[1][3], sigdigits=6), " seconds (on thread 1)")
    println("        finding redundant happairs     = ", round(quad_timers[1][4], sigdigits=6), " seconds (on thread 1)")
    println("    Phasing by dynamic programming  = ", round(phase_time, sigdigits=6), " seconds")
    println("    Imputation                      = ", round(impute_time, sigdigits=6), " seconds\n")

    return ph, redundant_haplotypes
end

"""
    phase!(X, H, width=400, verbose=true)

Phasing (haplotying) of genotype matrix `X` from a pool of haplotypes `H`
by dynamic programming. 

# Input
* `ph`: A vector of `HaplotypeMosaicPair` keeping track of each person's phase information.
* `X`: `p x n` matrix with missing values. Each column is genotypes of an individual.
* `compressed_Hunique`: A `CompressedHaplotypes` keeping track of unique haplotypes for each window and some other information
* `redundant_haplotypes`: Vector of optimal haplotype pairs across windows. The haplotype pairs are indices to the full haplotype set and NOT the compressed haplotypes
* `chunk_offset`: Shifts SNPs if a chromosome had been chunked. (not currently implemented)
"""
function phase!(
    ph::Vector{HaplotypeMosaicPair},
    X::AbstractMatrix{Union{Missing, T}},
    compressed_Hunique::CompressedHaplotypes,
    redundant_haplotypes::Vector{Vector{Vector{Tuple{Int, Int}}}},
    X_pos::Vector{Int};
    chunk_offset::Int = 0,
    ) where T <: Real

    # declare some constants
    people = size(X, 2)
    haplotypes = nhaplotypes(compressed_Hunique)
    snps = size(X, 1)
    width = compressed_Hunique.width
    windows = floor(Int, snps / width)
    H_pos = compressed_Hunique.pos
    last_window_width = snps - (windows - 1) * width 
    XtoH_idx = indexin(X_pos, H_pos)

    # allocate working arrays
    sol_path = [Vector{Tuple{Int, Int}}(undef, windows) for i in 1:Threads.nthreads()]
    nxt_pair = [[Int[] for i in 1:windows] for i in 1:Threads.nthreads()]
    tree_err = [[Float64[] for i in 1:windows] for i in 1:Threads.nthreads()]
    pmeter   = Progress(people, 5, "Merging breakpoints...")

    # loop over each person
    # first  1/3: ((w - 2) * width + 1):((w - 1) * width)
    # middle 1/3: ((w - 1) * width + 1):(      w * width)
    # last   1/3: (      w * width + 1):((w + 1) * width)
    for i in 1:people
        # first find optimal haplotype pair in each window using dynamic programming
        id = Threads.threadid()
        connect_happairs!(sol_path[id], nxt_pair[id], tree_err[id], redundant_haplotypes[i], λ = 1.0)

        # phase first window 
        h1 = complete_idx_to_unique_all_idx(sol_path[id][1][1], 1, compressed_Hunique)
        h2 = complete_idx_to_unique_all_idx(sol_path[id][1][2], 1, compressed_Hunique)
        push!(ph[i].strand1.start, 1 + chunk_offset)
        push!(ph[i].strand1.window, 1) 
        push!(ph[i].strand1.haplotypelabel, h1)
        push!(ph[i].strand2.start, 1 + chunk_offset)
        push!(ph[i].strand2.window, 1)
        push!(ph[i].strand2.haplotypelabel, h2)

        # don't search breakpoints
        for w in 2:windows
            u, j = sol_path[id][w - 1] # haplotype pair in previous window
            k, l = sol_path[id][w]     # haplotype pair in current window

            # switch current window's pair order if 1 or 2 haplotype match
            if (u == l && j == k) || (j == k && u ≠ l) || (u == l && j ≠ k)
                k, l = l, k 
                sol_path[id][w] = (k, l)
            end

            # map hap1 and hap2 back to unique index in given window
            h1 = complete_idx_to_unique_all_idx(k, w, compressed_Hunique)
            h2 = complete_idx_to_unique_all_idx(l, w, compressed_Hunique)

            push!(ph[i].strand1.start, chunk_offset + (w - 1) * width + 1)
            push!(ph[i].strand1.haplotypelabel, h1)
            push!(ph[i].strand1.window, w)
            push!(ph[i].strand2.start, chunk_offset + (w - 1) * width + 1)
            push!(ph[i].strand2.haplotypelabel, h2)
            push!(ph[i].strand2.window, w)
        end

        # search breakpoints 
        # for w in 2:windows
        #     # get genotype vector spanning 2 windows
        #     Xwi_start = (w - 2) * width + 1
        #     Xwi_end = (w == windows ? snps : w * width)
        #     Xwi = view(X, Xwi_start:Xwi_end, i)

        #     # find optimal breakpoint if there is one
        #     sol_path[id][w], bkpts = continue_haplotype(Xwi, compressed_Hunique, 
        #         w, sol_path[id][w - 1], sol_path[id][w])

        #     # record strand 1 info
        #     update_phase!(ph[i].strand1, compressed_Hunique, bkpts[1], sol_path[id][w - 1][1], 
        #         sol_path[id][w][1], w, width, chunk_offset, XtoH_idx, Xwi_start, Xwi_end)
        #     # record strand 2 info
        #     update_phase!(ph[i].strand2, compressed_Hunique, bkpts[2], sol_path[id][w - 1][2], 
        #         sol_path[id][w][2], w, width, chunk_offset, XtoH_idx, Xwi_start, Xwi_end)
        # end

        # update progress
        next!(pmeter)
    end
end

"""
Helper function for updating phase information after breakpoints have been identified
between windows `w - 1` and `w`. Every window have 0 or 1 breakpoint. Here indices in 
`ph.start` are recorded in terms of X's index.  

Caveat: technically it is possible for a window to have 2 breakpoints (that might even overlap)
since we can have the previous and next window both extend into the current one, but hopefully
this is extremely rare.
"""
function update_phase!(ph::HaplotypeMosaic, compressed_Hunique::CompressedHaplotypes,
    bkpt::Int, hap_prev::Int, hap_curr::Int, w::Int, width::Int, chunk_offset::Int,
    XtoH_idx::AbstractVector, Xwi_start::Int, Xwi_end::Int)

    # no breakpoints
    if bkpt == -1
        h = complete_idx_to_unique_all_idx(hap_curr, w, compressed_Hunique)
        push!(ph.start, chunk_offset + (w - 1) * width + 1)
        push!(ph.haplotypelabel, h)
        push!(ph.window, w)
        return nothing
    end

    # previous window's haplotype completely covers current window 
    if bkpt == length(Xwi_start:Xwi_end)
        h = complete_idx_to_unique_all_idx(hap_prev, w, compressed_Hunique)
        push!(ph.start, chunk_offset + (w - 1) * width + 1)
        push!(ph.haplotypelabel, h)
        push!(ph.window, w)
        return nothing
    end

    X_bkpt_end = Xwi_start + bkpt
    Xwi_mid = (w - 1) * width + 1

    if Xwi_mid <= X_bkpt_end <= Xwi_end
        # previous window extends to current window 
        h1 = complete_idx_to_unique_all_idx(hap_prev, w, compressed_Hunique)
        push!(ph.start, chunk_offset + Xwi_mid)
        push!(ph.haplotypelabel, h1)
        push!(ph.window, w)
        # 2nd part of current window
        h2 = complete_idx_to_unique_all_idx(hap_curr, w, compressed_Hunique)
        push!(ph.start, chunk_offset + X_bkpt_end)
        push!(ph.haplotypelabel, h2)
        push!(ph.window, w)
    elseif X_bkpt_end < Xwi_mid
        # current window extends to previous window
        h1 = complete_idx_to_unique_all_idx(hap_curr, w - 1, compressed_Hunique)
        push!(ph.start, chunk_offset + X_bkpt_end)
        push!(ph.haplotypelabel, h1)
        push!(ph.window, w - 1)
        # update current window
        h2 = complete_idx_to_unique_all_idx(hap_curr, w, compressed_Hunique)
        push!(ph.start, chunk_offset + Xwi_mid)
        push!(ph.haplotypelabel, h2)
        push!(ph.window, w)
    else
        # println("H_bkpt_pos = $H_bkpt_pos, Hw_start=$Hw_start, Hw_mid=$Hw_mid, Hw_end=$Hw_end ")
        error("update_phase!: bkpt does not satisfy -1 <= bkpt <= 2width! Shouldn't be possible")
    end

    return nothing
end
