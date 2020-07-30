"""
    phase(tgtfile, reffile; [outfile], [impute], [width], [recreen], 
    [thinning_factor], [dynamic_programming])

Main function of MendelImpute program. Phasing (haplotying) of `tgtfile` from a
pool of haplotypes `reffile` by sliding windows and saves result in `outfile`.

# Input
- `tgtfile`: VCF or PLINK files. VCF files should end in `.vcf` or `.vcf.gz`.
    PLINK files should exclude `.bim/.bed/.fam` suffixes but the trio must all be present in the directory.
- `reffile`: VCF or compressed Julia binary files. VCF files should end in 
    `.vcf` or `.vcf.gz`. Acceptable Julia binary formats includes `.jld2`
    (fastest read time) and `.jlso` (smallest file size).

# Optional Inputs
- `outfile`: output filename ending in `.vcf.gz` or `.vcf`. Output genotypes
    will have no missing data.
- `impute`: If `true`, untyped SNPs will be imputed, otherwise only missing snps
    in `tgtfile` will be imputed.  (default `false`)
- `phase`: If `true`, all output genotypes will be phased. Otherwise all output
    genotypes will be unphased.
- `width`: number of SNPs (markers) in each haplotype window. (default `512`)
- `rescreen`: This option saves a number of top haplotype pairs when solving the
    least squares objective, and re-minimize least squares on just observed data.
- `thinning_factor`: This option solves the least squares objective on only
    "thining_factor" unique haplotypes.
"""
function phase(
    tgtfile::AbstractString,
    reffile::AbstractString;
    outfile::AbstractString = "imputed." * tgtfile,
    impute::Bool = true,
    phase::Bool = false,
    width::Int = 512,
    rescreen::Bool = false,
    max_haplotypes::Int = 800,
    thinning_factor::Union{Nothing, Int} = nothing,
    scale_allelefreq::Bool = false,
    dynamic_programming::Bool = true,
    lasso::Union{Nothing, Int} = nothing
    )

    if dynamic_programming
        error("Currently dynamic programming routine is broken! Sorry!")
    end
    if !impute
        error("Currently one cannot impute only typed SNPs! Sorry!")
    end

    # import reference data
    println("Importing reference haplotype data..."); flush(stdout)
    import_data_start = time()
    if endswith(reffile, ".jld2")
        @load reffile compressed_Hunique
        width == compressed_Hunique.width || error("Specified width = $width" *
            " does not equal $(compressed_Hunique.width) = width in .jdl2 file")
    elseif endswith(reffile, ".jlso")
        loaded = JLSO.load(reffile)
        compressed_Hunique = loaded[:compressed_Hunique]
        width == compressed_Hunique.width || error("Specified width = $width" *
            " does not equal $(compressed_Hunique.width) = width in .jlso file")
    elseif endswith(reffile, ".vcf") || endswith(reffile, ".vcf.gz")
        # compress and filter VCF files for unique haplotypes in each window
        @info "VCF files detected: compressing reference file to .jlso format.."
        compressed_Hunique = compress_haplotypes(reffile, tgtfile,
            "compressed." * reffile, width)
    else
        error("Unrecognized reference file format: only VCF (ends in" * 
            " .vcf or .vcf.gz), `.jlso`, or `.jld2` files are acceptable.")
    end
    # import genotype data
    if endswith(tgtfile, ".vcf") || endswith(tgtfile, ".vcf.gz")
        X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt = 
            VCFTools.convert_gt(UInt8, tgtfile, trans=true, 
            save_snp_info=true, msg = "Importing genotype file...")
    elseif isplink(tgtfile)
        # PLINK files
        X_snpdata = SnpArrays.SnpData(tgtfile)
        X = convert(Matrix{UInt8}, X_snpdata.snparray')
        X_sampleID = X_snpdata.person_info[!, :iid]
        X_chr = X_snpdata.snp_info[!, :chromosome]
        X_pos = X_snpdata.snp_info[!, :position]
        X_ids = X_snpdata.snp_info[!, :snpid]
        X_ref = X_snpdata.snp_info[!, :allele1]
        X_alt = X_snpdata.snp_info[!, :allele2]
    else
        error("Unrecognized target file format: target file can only be VCF" *
            " files (ends in .vcf or .vcf.gz) or PLINK files (do not include" *
            " .bim/bed/fam and all three files must exist in 1 directory)")
    end
    import_data_time = time() - import_data_start

    # some constants and timers
    people = size(X, 2)
    tgt_snps = size(X, 1)
    ref_snps = length(compressed_Hunique.pos)
    windows = floor(Int, tgt_snps / width)
    num_unique_haps = round(Int, avg_haplotypes_per_window(compressed_Hunique))

    # working arrays
    ph = [HaplotypeMosaicPair(ref_snps) for i in 1:people]
    haplotype1 = [zeros(Int32, windows) for i in 1:people]
    haplotype2 = [zeros(Int32, windows) for i in 1:people]
    # if dynamic_programming
    #     redundant_haplotypes = [[Tuple{Int32, Int32}[] for i in
    #         1:num_windows_per_chunks] for j in 1:people]
    #     [[sizehint!(redundant_haplotypes[j][i], 1000) for i in
    #         1:num_windows_per_chunks] for j in 1:people]
    # end

    #
    # find best happairs for each window
    #
    calculate_happairs_start = time()
    haptimers = compute_optimal_haplotypes!(haplotype1, haplotype2, 
        compressed_Hunique, X, X_pos, lasso, thinning_factor, scale_allelefreq, 
        max_haplotypes, rescreen)
    screen_flanking_windows!(haplotype1, haplotype2, compressed_Hunique, X)
    calculate_happairs_time = time() - calculate_happairs_start

    #
    # phasing (haplotyping) + breakpoint search
    #
    phase_start = time()
    if dynamic_programming
        phase!(ph, X, compressed_Hunique, redundant_haplotypes, X_pos,
        1:windows) # dynamic programming
    else
        phase_fast!(ph, X, compressed_Hunique, haplotype1, haplotype2, X_pos,
            1:windows) # phase window-by-window
    end
    phase_time = time() - phase_start

    #
    # impute step
    #
    impute_start = time()
    XtoH_idx = indexin(X_pos, compressed_Hunique.pos)
    if impute # imputes typed and untyped SNPs
        # convert phase's starting position from X's index to H's index
        update_marker_position!(ph, XtoH_idx)

        if phase # output genotypes all phased
            X1 = BitArray(undef, ref_snps, people)
            X2 = BitArray(undef, ref_snps, people)

            # impute and write to file
            impute!(X1, X2, compressed_Hunique, ph)
            write(outfile, X1, X2, compressed_Hunique, X_sampleID)
        else # output genotypes all unphased
            X_full = Matrix{Union{Missing, UInt8}}(missing, ref_snps, people)
            copyto!(@view(X_full[XtoH_idx, :]), X) # keep known entries

            # impute and write to file
            impute_discard_phase!(X_full, compressed_Hunique, ph)
            write(outfile, X_full, compressed_Hunique, X_sampleID)
        end
    else # impute only missing entries in typed SNPs
        if phase
            X1 = BitArray(undef, size(X, 1), size(X, 2))
            X2 = BitArray(undef, size(X, 1), size(X, 2))

            # impute and write to file
            impute!(X1, X2, compressed_Hunique, ph)
            write(outfile, X1, X2, compressed_Hunique, X_sampleID, XtoH_idx)
        else
            impute_discard_phase!(X, compressed_Hunique, ph)
            write(outfile, X, compressed_Hunique, X_sampleID, XtoH_idx)
        end
    end
    impute_time = time() - impute_start

    # print timing results
    println("Total windows = $windows, averaging ~ $num_unique_haps " *
        "unique haplotypes per window.\n")
    println("Timings: ")
    println("    Data import                     = ", 
        round(import_data_time, sigdigits=6), " seconds")
    println("    Computing haplotype pair        = ", 
        round(calculate_happairs_time, sigdigits=6), " seconds")
    haptimers[1] != 0 && println("        screening for top haplotypes   = ", 
        round(haptimers[1], sigdigits=6), " seconds per thread")
    println("        BLAS3 mul! to get M and N      = ", 
        round(haptimers[2], sigdigits=6), " seconds per thread")
    println("        haplopair search               = ", 
        round(haptimers[3], sigdigits=6), " seconds per thread")
    haptimers[4] != 0 && println("        min least sq on observed data  = ", 
        round(haptimers[4], sigdigits=6), " seconds per thread")
    println("        index conversion               = ", 
        round(haptimers[5], sigdigits=6), " seconds per thread")
    dynamic_programming ? println("    Phasing by dynamic programming  = ", 
                          round(phase_time, sigdigits=6), " seconds") :
                          println("    Phasing by win-win intersection = ", 
                          round(phase_time, sigdigits=6), " seconds")
    println("    Imputation                      = ", 
        round(impute_time, sigdigits=6), " seconds\n")

    return ph
end

isplink(tgtfile::AbstractString) = isfile(tgtfile * ".bed") && 
                                   isfile(tgtfile * ".fam") && 
                                   isfile(tgtfile * ".bim")

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
    redundant_haplotypes::Vector{Vector{Vector{Tuple{Int32, Int32}}}},
    X_pos::Vector{Int},
    total_window::Int,
    winrange::UnitRange
    ) where T <: Real

    # declare some constants
    people = size(X, 2)
    snps = size(X, 1)
    width = compressed_Hunique.width
    windows = length(winrange)
    H_pos = compressed_Hunique.pos
    chunk_offset = (first(winrange) - 1) * width

    # allocate working arrays
    sol_path = [Vector{Tuple{Int32, Int32}}(undef, windows) for i in 1:Threads.nthreads()]
    nxt_pair = [[Int32[] for i in 1:windows] for i in 1:Threads.nthreads()]
    tree_err = [[Float64[] for i in 1:windows] for i in 1:Threads.nthreads()]

    # loop over each person
    # first  1/3: ((w - 2) * width + 1):((w - 1) * width)
    # middle 1/3: ((w - 1) * width + 1):(      w * width)
    # last   1/3: (      w * width + 1):((w + 1) * width)
    ThreadPools.@qthreads for i in 1:people
        # first find optimal haplotype pair in each window using dynamic programming
        id = Threads.threadid()
        connect_happairs!(sol_path[id], nxt_pair[id], tree_err[id], redundant_haplotypes[i], λ = 1.0)

        # phase first window
        h1 = complete_idx_to_unique_all_idx(sol_path[id][1][1], first(winrange), compressed_Hunique)
        h2 = complete_idx_to_unique_all_idx(sol_path[id][1][2], first(winrange), compressed_Hunique)
        push!(ph[i].strand1.start, 1 + chunk_offset)
        push!(ph[i].strand1.window, first(winrange))
        push!(ph[i].strand1.haplotypelabel, h1)
        push!(ph[i].strand2.start, 1 + chunk_offset)
        push!(ph[i].strand2.window, first(winrange))
        push!(ph[i].strand2.haplotypelabel, h2)

        # don't search breakpoints
        # for w in 2:windows
        #     u, j = sol_path[id][w - 1] # haplotype pair in previous window
        #     k, l = sol_path[id][w]     # haplotype pair in current window

        #     # switch current window's pair order if 1 or 2 haplotype match
        #     if (u == l && j == k) || (j == k && u ≠ l) || (u == l && j ≠ k)
        #         k, l = l, k
        #         sol_path[id][w] = (k, l)
        #     end

        #     # map hap1 and hap2 back to unique index in given window
        #     h1 = complete_idx_to_unique_all_idx(k, w, compressed_Hunique)
        #     h2 = complete_idx_to_unique_all_idx(l, w, compressed_Hunique)

        #     push!(ph[i].strand1.start, chunk_offset + (w - 1) * width + 1)
        #     push!(ph[i].strand1.haplotypelabel, h1)
        #     push!(ph[i].strand1.window, w)
        #     push!(ph[i].strand2.start, chunk_offset + (w - 1) * width + 1)
        #     push!(ph[i].strand2.haplotypelabel, h2)
        #     push!(ph[i].strand2.window, w)
        # end

        # search breakpoints
        for (w, absolute_w) in zip(2:windows, winrange[2:windows])
            # get genotype vector spanning 2 windows
            Xwi_start = (absolute_w - 2) * width + 1
            Xwi_end = (absolute_w == total_window ? snps : absolute_w * width)
            Xwi = view(X, Xwi_start:Xwi_end, i)

            # find optimal breakpoint if there is one
            sol_path[id][w], bkpts = continue_haplotype(Xwi, compressed_Hunique,
                absolute_w, sol_path[id][w - 1], sol_path[id][w])

            # record strand 1 info
            update_phase!(ph[i].strand1, compressed_Hunique, bkpts[1], sol_path[id][w - 1][1],
                sol_path[id][w][1], absolute_w, width, Xwi_start, Xwi_end)
            # record strand 2 info
            update_phase!(ph[i].strand2, compressed_Hunique, bkpts[2], sol_path[id][w - 1][2],
                sol_path[id][w][2], absolute_w, width, Xwi_start, Xwi_end)
        end
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
    bkpt::Int, hap_prev, hap_curr, window::Int, width::Int, Xwi_start::Int, Xwi_end::Int)

    # no breakpoints
    if bkpt == -1
        h = complete_idx_to_unique_all_idx(hap_curr, window, compressed_Hunique)
        push!(ph.start, (window - 1) * width + 1)
        push!(ph.haplotypelabel, h)
        push!(ph.window, window)
        return nothing
    end

    # previous window's haplotype completely covers current window
    if bkpt == length(Xwi_start:Xwi_end)
        h = complete_idx_to_unique_all_idx(hap_prev, window, compressed_Hunique)
        push!(ph.start, (window - 1) * width + 1)
        push!(ph.haplotypelabel, h)
        push!(ph.window, window)
        return nothing
    end

    X_bkpt_end = Xwi_start + bkpt
    Xwi_mid = (window - 1) * width + 1

    if Xwi_mid <= X_bkpt_end <= Xwi_end
        # previous window extends to current window
        h1 = complete_idx_to_unique_all_idx(hap_prev, window, compressed_Hunique)
        push!(ph.start, Xwi_mid)
        push!(ph.haplotypelabel, h1)
        push!(ph.window, window)
        # 2nd part of current window
        h2 = complete_idx_to_unique_all_idx(hap_curr, window, compressed_Hunique)
        push!(ph.start, X_bkpt_end)
        push!(ph.haplotypelabel, h2)
        push!(ph.window, window)
    elseif X_bkpt_end < Xwi_mid
        # current window extends to previous window
        h1 = complete_idx_to_unique_all_idx(hap_curr, window - 1, compressed_Hunique)
        push!(ph.start, X_bkpt_end)
        push!(ph.haplotypelabel, h1)
        push!(ph.window, window - 1)
        # update current window
        h2 = complete_idx_to_unique_all_idx(hap_curr, window, compressed_Hunique)
        push!(ph.start, Xwi_mid)
        push!(ph.haplotypelabel, h2)
        push!(ph.window, window)
    else
        # println("H_bkpt_pos = $H_bkpt_pos, Hw_start=$Hw_start, Hw_mid=$Hw_mid, Hw_end=$Hw_end ")
        error("update_phase!: bkpt does not satisfy -1 <= bkpt <= 2width! Shouldn't be possible")
    end

    return nothing
end

function phase_fast!(
    ph::Vector{HaplotypeMosaicPair},
    X::AbstractMatrix{Union{Missing, T}},
    compressed_Hunique::CompressedHaplotypes,
    haplotype1::AbstractVector,
    haplotype2::AbstractVector,
    X_pos::Vector{Int},
    winrange::UnitRange
    ) where T <: Real

    # declare some constants
    people = size(X, 2)
    snps = size(X, 1)
    haplotypes = nhaplotypes(compressed_Hunique)
    width = compressed_Hunique.width
    windows = length(winrange)
    H_pos = compressed_Hunique.pos

    # working arrays
    seen = [BitSet() for _ in 1:Threads.nthreads()]
    survivors1 = [Int32[] for _ in 1:Threads.nthreads()]
    survivors2 = [Int32[] for _ in 1:Threads.nthreads()]
    for id in 1:Threads.nthreads()
        sizehint!(seen[id], haplotypes)
        sizehint!(survivors1[id], haplotypes)
        sizehint!(survivors2[id], haplotypes)
    end

    # first  1/3: ((w - 2) * width + 1):((w - 1) * width)
    # middle 1/3: ((w - 1) * width + 1):(      w * width)
    # last   1/3: (      w * width + 1):((w + 1) * width)
    ThreadPools.@qthreads for i in 1:people
        id = Threads.threadid()

        # First pass to phase each sample window-by-window
        phase_sample!(haplotype1[i], haplotype2[i], compressed_Hunique,
            seen[id], survivors1[id], survivors2[id])

        # record info for first window
        hap1 = haplotype1[i][1] # complete idx
        hap2 = haplotype2[i][1] # complete idx
        h1 = complete_idx_to_unique_all_idx(hap1, 1,
            compressed_Hunique) #unique haplotype idx in window 1
        h2 = complete_idx_to_unique_all_idx(hap2, 1,
            compressed_Hunique) #unique haplotype idx in window 1
        push!(ph[i].strand1.start, 1)
        push!(ph[i].strand1.window, 1)
        push!(ph[i].strand1.haplotypelabel, h1)
        push!(ph[i].strand2.start, 1)
        push!(ph[i].strand2.window, 1)
        push!(ph[i].strand2.haplotypelabel, h2)

        # Second pass to find optimal break points and record info to phase
        @inbounds for w in 2:windows
            # get genotype vector spanning 2 windows
            Xwi_start = (w - 2) * width + 1
            Xwi_end = (w == windows ? snps : w * width)
            Xwi = view(X, Xwi_start:Xwi_end, i)

            # previous and current haplotypes for both strands
            hap1_prev = haplotype1[i][w - 1]
            hap2_prev = haplotype2[i][w - 1]
            hap1_curr = haplotype1[i][w]
            hap2_curr = haplotype2[i][w]

            # find optimal breakpoint if there is one
            _, bkpts = continue_haplotype(Xwi, compressed_Hunique,
                w, (hap1_prev, hap2_prev), (hap1_curr, hap2_curr))

            # record strand 1 info
            update_phase!(ph[i].strand1, compressed_Hunique, bkpts[1],
                hap1_prev, hap1_curr, w, width, Xwi_start, Xwi_end)
            # record strand 2 info
            update_phase!(ph[i].strand2, compressed_Hunique, bkpts[2],
                hap2_prev, hap2_curr, w, width, Xwi_start, Xwi_end)
        end
    end
end
