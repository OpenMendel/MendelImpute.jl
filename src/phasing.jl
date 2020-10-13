###### This file is part of the MendelImpute.jl package.
###### This is the main entry point of the MendelImpute package. 
###### It also contains code to phase window-by-window. 

"""
    phase(tgtfile::String, reffile::String, outfile::String; [impute::Bool],
        [phase::Bool], [dosage::Bool], [recreen::Bool], [max_haplotypes::Int], 
        [stepwise::Int], [thinning_factor::Int], [scale_allelefreq::Bool], 
        [dynamic_programming::Bool])

Main function of MendelImpute program. Phasing (haplotying) of `tgtfile` from a
pool of haplotypes `reffile` by sliding windows and saves result in `outfile`.
All SNPs in `tgtfile` must be present in `reffile`. Per-sample imputation score
(lower is better) will be saved in a file ending in `sample.error`.

# Input
- `tgtfile`: VCF or PLINK files. VCF files should end in `.vcf` or `.vcf.gz`.
    PLINK files should exclude `.bim/.bed/.fam` suffixes but the trio must all
    be present in the same directory.
- `reffile`: Reference haplotype file ending in `.vcf`, `.vcf.gz`, or `.jlso` 
    (compressed binary files).
- `outfile`: output filename ending in `.vcf.gz`, `.vcf`, or `.jlso`. VCF output
    genotypes will have no missing data. If ending in `.jlso`, will output
    ultra-compressed data structure recording `HaplotypeMosaicPair`s for 
    each sample.

# Optional Inputs
- `impute`: If `true`, imputes every SNPs in `reffile` to `tgtfile`. Otherwise
    only missing snps in `tgtfile` will be imputed.
- `phase`: If `true`, all output genotypes will be phased, but observed data 
    (minor allele count) may be changed. If `phase=false` all output genotypes
    will be unphased but observed minor allele count will not change.
- `dosage`: If `true`, will assume target matrix are dosages for imputation. Note
    this means the genotype matrix will be entirely single precision. 
- `rescreen`: This option is more computationally intensive but gives more
    accurate results. It saves a number of top haplotype pairs when solving
    the least squares objective, and re-minimize least squares on just
    observed data.
- `max_haplotypes`: Maximum number of haplotypes for using global search. Windows
    exceeding this number of unique haplotypes will be searched using a
    heuristic. A non-zero `stepscreen` or `thinning_factor` need to be specified 
- `stepwise`: If an integer is specified, will solve the least squares objective
    by first finding `stepwise` top haplotypes using a stepwise heuristic then
    finds the next haplotype using global search. Uses `max_haplotypes`. 
- `thinning_factor`: If an integer is specified, will solve the least squares
    objective on only `thining_factor` unique haplotypes. Uses `max_haplotypes`.
- `scale_allelefreq`: Boolean indicating whether to give rare SNPs more weight
    scaled by `wᵢ = 1 / √2p(1-p)` where max weight is 2. 
- `dynamic_programming`: Boolean indicating whether to phase with a global 
    search that finds the longest haplotype stretch over all windows. (Currently
    broken, sorry!)
"""
function phase(
    tgtfile::AbstractString,
    reffile::AbstractString,
    outfile::AbstractString;
    impute::Bool = true,
    phase::Bool = true,
    dosage::Bool = false,
    rescreen::Bool = false,
    max_haplotypes::Int = 800,
    stepwise::Union{Nothing, Int} = nothing,
    thinning_factor::Union{Nothing, Int} = nothing,
    scale_allelefreq::Bool = false,
    dynamic_programming::Bool = false
    )

    # first handle errors
    if dynamic_programming
        error("Currently dynamic programming routine is broken! Sorry!")
    end
    endswith(outfile, ".jlso") || endswith(outfile, ".vcf") || 
        endswith(outfile, ".vcf.gz") || error("Output file name must end with" * 
        " .vcf or .vcf.gz or .jlso!")
    ultra_compress = endswith(outfile, ".jlso") ? true : false
    ultra_compress && !phase && error("Ultra compressed output must be phased!")

    # import reference data
    println("Number of threads = ", Threads.nthreads())
    println("Importing reference haplotype data..."); flush(stdout)
    ref_import_start = time()
    if endswith(reffile, ".jlso")
        loaded = JLSO.load(reffile)
        compressed_Hunique = loaded[:compressed_Hunique]
    elseif endswith(reffile, ".vcf") || endswith(reffile, ".vcf.gz")
        error("reference panel is in VCF for: please first compress reference" *
            " file to .jlso format using the compress_haplotypes() function.")
    elseif endswith(reffile, ".jld2")
        @load reffile compressed_Hunique 
    else
        error("Unrecognized reference file format: only VCF (ends in" * 
            " .vcf or .vcf.gz) or `.jlso` files are acceptable.")
    end
    ref_import_time = time() - ref_import_start

    # import genotype data
    genotype_import_start = time()
    if (endswith(tgtfile, ".vcf") || endswith(tgtfile, ".vcf.gz")) && !dosage
        X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt = 
            VCFTools.convert_gt(UInt8, tgtfile, trans=true, 
            save_snp_info=true, msg = "Importing genotype file...")
    elseif (endswith(tgtfile, ".vcf") || endswith(tgtfile, ".vcf.gz")) && dosage
        X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt = 
            VCFTools.convert_ds(Float32, tgtfile, trans=true, 
            save_snp_info=true, msg = "Importing genotype file as dosages...")
    elseif isplink(tgtfile)
        dosage && error("PLINK files detected but dosage = true!")
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
    genotype_import_time = time() - genotype_import_start
    import_data_time = time() - ref_import_start

    # some constants and timers
    people = size(X, 2)
    tgt_snps = size(X, 1)
    ref_snps = length(compressed_Hunique.pos)
    windows = nwindows(compressed_Hunique)
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
        compressed_Hunique, X, X_pos, stepwise, thinning_factor,
        scale_allelefreq, max_haplotypes, rescreen)
    # screen_flanking_windows!(haplotype1, haplotype2, compressed_Hunique, X)
    calculate_happairs_time = time() - calculate_happairs_start

    #
    # phasing (haplotyping) + breakpoint search
    #
    phase_start = time()
    if dynamic_programming
        phase!(ph, X, compressed_Hunique, redundant_haplotypes, X_pos,
            1:windows)
    elseif ultra_compress # phase window-by-window for outputing ultra compressed format
        phasetimers = phase_fast_compressed!(ph, X, compressed_Hunique, 
            haplotype1, haplotype2)
    else # phase window-by-window for outputing VCF files
        phasetimers = phase_fast!(ph, X, compressed_Hunique, haplotype1,
            haplotype2, impute)
    end
    phase_time = time() - phase_start

    #
    # impute step
    #
    impute_start = time()
    write_time = 0.0
    XtoH_idx = indexin(X_pos, compressed_Hunique.pos)
    # get each snp's imputation score
    snpscore = typed_snpscore(X, ph, compressed_Hunique)
    if impute # imputes typed and untyped SNPs
        # get quality score for untyped SNPs
        complete_snpscore = untyped_snpscore(ref_snps, snpscore, XtoH_idx)

        # convert phase's starting position from X's index to H's index
        update_marker_position!(ph, XtoH_idx)

        if ultra_compress # output ultra-compressed, phased genotypes in 
            write_time += @elapsed JLSO.save(outfile, :ph => ph, 
                :sampleID => X_sampleID, format=:julia_serialize, 
                compression=:gzip)
        elseif phase # output genotypes all phased
            X1 = BitArray(undef, ref_snps, people)
            X2 = BitArray(undef, ref_snps, people)

            # impute and write to file
            impute!(X1, X2, compressed_Hunique, ph, impute)
            write_time += @elapsed write(outfile, (X1, X2), compressed_Hunique, 
                X_sampleID, complete_snpscore, XtoH_idx)
        else # output genotypes all unphased
            X_full = Matrix{Union{Missing, UInt8}}(missing, ref_snps, people)
            copyto!(@view(X_full[XtoH_idx, :]), X) # keep known entries

            # impute and write to file
            impute_discard_phase!(X_full, compressed_Hunique, ph, impute)
            write_time += @elapsed write(outfile, X_full, compressed_Hunique, 
                X_sampleID, complete_snpscore, XtoH_idx)
        end
    else # impute only missing entries in typed SNPs
        if ultra_compress # output ultra-compressed, phased genotypes in 
            write_time += @elapsed JLSO.save(outfile, :ph => ph, 
                :sampleID => X_sampleID, format=:julia_serialize, 
                compression=:gzip)
        elseif phase
            X1 = BitArray(undef, size(X, 1), size(X, 2))
            X2 = BitArray(undef, size(X, 1), size(X, 2))

            # impute and write to file
            impute!(X1, X2, compressed_Hunique, ph, impute)
            write_time += @elapsed write(outfile, (X1, X2), compressed_Hunique, 
                X_sampleID, snpscore, XtoH_idx, false)
        else
            impute_discard_phase!(X, compressed_Hunique, ph, impute)
            write_time += @elapsed write(outfile, X, compressed_Hunique, 
                X_sampleID, snpscore, XtoH_idx, false)
        end
    end
    # write per-sample error to output file
    if endswith(outfile, ".jlso")
        strip_chr = 4
    elseif endswith(outfile, ".vcf")
        strip_chr = 3
    else # .vcf.gz
        strip_chr = 6
    end
    # error_filename = outfile[1:end-strip_chr]
    # write_time += @elapsed begin
    #     open(error_filename * "sample.error", "w") do io
    #         print(io, "ID,error\n")
    #         for i in eachindex(X_sampleID)
    #             @inbounds print(io, X_sampleID[i], ",", sum(haploscore[i]), "\n")
    #         end
    #     end
    # end
    impute_time = time() - impute_start
    impute_nonwrite_time = impute_time - write_time

    #
    # print timing results
    #
    println("Total windows = $windows, averaging ~ $num_unique_haps " *
        "unique haplotypes per window.\n")
    println("Timings: ")
    println("    Data import                     = ", 
        round(import_data_time, sigdigits=6), " seconds")
    println("        import target data             = ", 
        round(genotype_import_time, sigdigits=6), " seconds")
    println("        import compressed haplotypes   = ", 
        round(ref_import_time, sigdigits=6), " seconds")
    println("    Computing haplotype pair        = ", 
        round(calculate_happairs_time, sigdigits=6), " seconds")
    haptimers[1*8] != 0 && println("        screening for top haplotypes   = ", 
        round(haptimers[1*8], sigdigits=6), " seconds per thread")
    println("        BLAS3 mul! to get M and N      = ", 
        round(haptimers[2*8], sigdigits=6), " seconds per thread")
    println("        haplopair search               = ", 
        round(haptimers[3*8], sigdigits=6), " seconds per thread")
    haptimers[4*8] != 0 && println("        min least sq on observed data  = ", 
        round(haptimers[4*8], sigdigits=6), " seconds per thread")
    println("        initializing missing           = ", 
        round(haptimers[5*8], sigdigits=6), " seconds per thread")
    println("        allocating and viewing         = ", 
        round(haptimers[6*8], sigdigits=6), " seconds per thread")
    println("        index conversion               = ", 
        round(haptimers[7*8], sigdigits=6), " seconds per thread")
    dynamic_programming ? println("    Phasing by dynamic programming  = ", 
                          round(phase_time, sigdigits=6), " seconds") :
                          println("    Phasing by win-win intersection = ", 
                          round(phase_time, sigdigits=6), " seconds")
    if !dynamic_programming
        println("        Window-by-window intersection  = ", 
            round(phasetimers[1*8], sigdigits=6), " seconds per thread")
        println("        Breakpoint search              = ", 
            round(phasetimers[2*8], sigdigits=6), " seconds per thread")
        println("        Recording result               = ", 
            round(phasetimers[3*8], sigdigits=6), " seconds per thread")
    end
    println("    Imputation                     = ", 
        round(impute_time, sigdigits=6), " seconds")
    println("        Imputing missing               = ", 
        round(impute_nonwrite_time, sigdigits=6), " seconds")
    println("        Writing to file                = ", 
        round(write_time, sigdigits=6), " seconds\n")

    total_time = time() - ref_import_start
    println("    Total time                      = ", 
        round(total_time, sigdigits=6), " seconds\n")

    return ph
end

isplink(tgtfile::AbstractString) = isfile(tgtfile * ".bed") && 
                                   isfile(tgtfile * ".fam") && 
                                   isfile(tgtfile * ".bim")

# """
#     phase!(X, H, width=400, verbose=true)

# Phasing (haplotying) of genotype matrix `X` from a pool of haplotypes `H`
# by dynamic programming.

# # Input
# * `ph`: A vector of `HaplotypeMosaicPair` keeping track of each person's
#     phaseinformation.
# * `X`: `p x n` matrix with missing values. Each column is genotypes of
#     an individual.
# * `compressed_Hunique`: A `CompressedHaplotypes` keeping track of unique
#     haplotypes for each window and some other information
# * `redundant_haplotypes`: Vector of optimal haplotype pairs across windows.
#     The haplotype pairs are indices to the full haplotype set and NOT the
#     compressed haplotypes
# """
# function phase!(
#     ph::Vector{HaplotypeMosaicPair},
#     X::AbstractMatrix{Union{Missing, T}},
#     compressed_Hunique::CompressedHaplotypes,
#     redundant_haplotypes::Vector{Vector{Vector{Tuple{Int32, Int32}}}},
#     X_pos::Vector{Int},
#     total_window::Int,
#     winrange::UnitRange
#     ) where T <: Real

#     # declare some constants
#     people = size(X, 2)
#     snps = size(X, 1)
#     width = compressed_Hunique.width
#     windows = length(winrange)
#     chunk_offset = (first(winrange) - 1) * width

#     # allocate working arrays
#     sol_path = [Vector{Tuple{Int32, Int32}}(undef, windows) for i in
#         1:Threads.nthreads()]
#     nxt_pair = [[Int32[] for i in 1:windows] for i in 1:Threads.nthreads()]
#     tree_err = [[Float64[] for i in 1:windows] for i in 1:Threads.nthreads()]

#     # loop over each person
#     # first  1/3: ((w - 2) * width + 1):((w - 1) * width)
#     # middle 1/3: ((w - 1) * width + 1):(      w * width)
#     # last   1/3: (      w * width + 1):((w + 1) * width)
#     ThreadPools.@qthreads for i in 1:people
#         # first find optimal haplotype pair in each window
#         id = Threads.threadid()
#         connect_happairs!(sol_path[id], nxt_pair[id], tree_err[id],
#             redundant_haplotypes[i], λ = 1.0)

#         # phase first window
#         h1 = complete_idx_to_unique_all_idx(sol_path[id][1][1],
#             first(winrange), compressed_Hunique)
#         h2 = complete_idx_to_unique_all_idx(sol_path[id][1][2],
#             first(winrange), compressed_Hunique)
#         push!(ph[i].strand1.start, 1 + chunk_offset)
#         push!(ph[i].strand1.window, first(winrange))
#         push!(ph[i].strand1.haplotypelabel, h1)
#         push!(ph[i].strand2.start, 1 + chunk_offset)
#         push!(ph[i].strand2.window, first(winrange))
#         push!(ph[i].strand2.haplotypelabel, h2)

#         # don't search breakpoints
#         # for w in 2:windows
#         #     u, j = sol_path[id][w - 1] # haplotype pair in previous window
#         #     k, l = sol_path[id][w]     # haplotype pair in current window

#         #     # switch current window's pair order if 1 or 2 haplotype match
#         #     if (u == l && j == k) || (j == k && u ≠ l) || (u == l && j ≠ k)
#         #         k, l = l, k
#         #         sol_path[id][w] = (k, l)
#         #     end

#         #     # map hap1 and hap2 back to unique index in given window
#         #     h1 = complete_idx_to_unique_all_idx(k, w, compressed_Hunique)
#         #     h2 = complete_idx_to_unique_all_idx(l, w, compressed_Hunique)

#         #     push!(ph[i].strand1.start, chunk_offset + (w - 1) * width + 1)
#         #     push!(ph[i].strand1.haplotypelabel, h1)
#         #     push!(ph[i].strand1.window, w)
#         #     push!(ph[i].strand2.start, chunk_offset + (w - 1) * width + 1)
#         #     push!(ph[i].strand2.haplotypelabel, h2)
#         #     push!(ph[i].strand2.window, w)
#         # end

#         # search breakpoints
#         for (w, absolute_w) in zip(2:windows, winrange[2:windows])
#             # get genotype vector spanning 2 windows
#             Xwi_start = (absolute_w - 2) * width + 1
#             Xwi_end = (absolute_w == total_window ? snps : absolute_w * width)
#             Xwi = view(X, Xwi_start:Xwi_end, i)

#             # find optimal breakpoint if there is one
#             sol_path[id][w], bkpts = continue_haplotype(Xwi, compressed_Hunique,
#                 absolute_w, sol_path[id][w - 1], sol_path[id][w])

#             # record strand 1 info
#             update_phase!(ph[i].strand1, compressed_Hunique, bkpts[1],
#                 sol_path[id][w - 1][1], sol_path[id][w][1], absolute_w, width,
#                 Xwi_start, Xwi_end)
#             # record strand 2 info
#             update_phase!(ph[i].strand2, compressed_Hunique, bkpts[2],
#                 sol_path[id][w - 1][2], sol_path[id][w][2], absolute_w, width,
#                 Xwi_start, Xwi_end)
#         end
#     end
# end

"""
    update_phase!(ph, compressed_Hunique, bkpt, ...)

Updates a person's phase information after breakpoints have been
identified between windows `w - 1` and `w`. Every window have 0 or 1
breakpoint. Here indices in `ph.start` are recorded in terms of X's index.
"""
function update_phase!(ph::HaplotypeMosaic,
    compressed_Hunique::CompressedHaplotypes, bkpt::Int64, hap_prev::Int32, 
    hap_curr::Int32, window::Int32, Xwi_start::Int64, Xwi_mid::Int64, 
    Xwi_end::Int64, impute_untyped::Bool)

    # specify conversion function
    convert = impute_untyped ? complete_idx_to_unique_all_idx :
        complete_idx_to_unique_typed_idx

    # no breakpoints or not searching double breakpoints
    if bkpt == -1 || bkpt == -2
        h = convert(hap_curr, window, compressed_Hunique)
        push_Mosaic!(ph, (Xwi_mid, h, window))
        return nothing
    end

    # previous window's haplotype completely covers current window
    if bkpt == length(Xwi_start:Xwi_end)
        h = convert(hap_prev, window, compressed_Hunique)
        push_Mosaic!(ph, (Xwi_mid, h, window))
        return nothing
    end

    X_bkpt_end = Xwi_start + bkpt

    if Xwi_mid <= X_bkpt_end <= Xwi_end
        # previous window extends to current window
        h1 = convert(hap_prev, window, compressed_Hunique)
        push_Mosaic!(ph, (Xwi_mid, h1, window))
        # 2nd part of current window
        h2 = convert(hap_curr, window, compressed_Hunique)
        push_Mosaic!(ph, (X_bkpt_end, h2, window))
    elseif X_bkpt_end < Xwi_mid
        # current window extends to previous window
        h1 = convert(hap_curr, Int32(window - 1), compressed_Hunique)
        push_Mosaic!(ph, (X_bkpt_end, h1, Int32(window - 1)))
        # update current window
        h2 = convert(hap_curr, window, compressed_Hunique)
        push_Mosaic!(ph, (Xwi_mid, h2, window))
    else
        error("update_phase!: bkpt does not satisfy -1 <= bkpt <= 2width!")
    end

    return nothing
end

"""
    phase_fast!(ph, X, compressed_Hunique, haplotype1, haplotype2, [impute_untyped])

Given optimal haplotype pairs in each window, performs window-by-window 
intersection heuristic to find longest spanning haplotypes (phasing) then
searches for optimal breakpoint. 

# Arguments
* `ph`: A vector of `HaplotypeMosaicPair` keeping track of each person's
    phaseinformation.
* `X`: `p x n` matrix with missing values. Each column is genotypes of
    an individual.
* `compressed_Hunique`: A `CompressedHaplotypes` keeping track of unique
    haplotypes for each window and some other information
* `haplotype1`: `haplotype1[w]` stores a optimal haplotype for window `w`. 
* `haplotype2`: `haplotype2[w]` stores a optimal haplotype for window `w`. 
* `impute_untyped`: Bool indicating whether untyped SNPs should be imputed. 

# Timers:
- `t1` = Window-by-window intersection
- `t2` = Breakpoint search
- `t3` = Recording result
"""
function phase_fast!(
    ph::Vector{HaplotypeMosaicPair},
    X::AbstractMatrix{Union{Missing, T}},
    compressed_Hunique::CompressedHaplotypes,
    haplotype1::AbstractVector,
    haplotype2::AbstractVector,
    impute_untyped::Bool
    ) where T <: Real

    # declare some constants
    people = size(X, 2)
    snps = size(X, 1)
    haplotypes = nhaplotypes(compressed_Hunique)
    winranges = compressed_Hunique.X_window_range
    windows = length(haplotype1[1])

    # working arrays
    survivors1 = [Int32[] for _ in 1:Threads.nthreads()]
    survivors2 = [Int32[] for _ in 1:Threads.nthreads()]
    for id in 1:Threads.nthreads()
        sizehint!(survivors1[id], haplotypes)
        sizehint!(survivors2[id], haplotypes)
    end
    timers = [zeros(3*8) for _ in 1:Threads.nthreads()] # 8 for spacing
    pmeter = Progress(people, 5, "Phasing...")

    # phase person by person
    # for i in 1:people
    Threads.@threads for i in 1:people
        id = Threads.threadid()

        # First pass to phase each sample window-by-window
        timers[id][8] += @elapsed phase_sample!(haplotype1[i], haplotype2[i],
            compressed_Hunique, survivors1[id], survivors2[id])

        # record info for first window
        timers[id][24] += @elapsed begin
            hap1 = haplotype1[i][1] # complete idx
            hap2 = haplotype2[i][1] # complete idx
            h1 = impute_untyped ? complete_idx_to_unique_all_idx(hap1, 1, 
                compressed_Hunique) : 
                complete_idx_to_unique_typed_idx(hap1, 1, compressed_Hunique)
            h2 = impute_untyped ? complete_idx_to_unique_all_idx(hap2, 1, 
                compressed_Hunique) : 
                complete_idx_to_unique_typed_idx(hap2, 1, compressed_Hunique)
            push!(ph[i].strand1.start, 1)
            push!(ph[i].strand1.haplotypelabel, h1)
            push!(ph[i].strand2.start, 1)
            push!(ph[i].strand2.haplotypelabel, h2)
            push!(ph[i].strand1.window, 1)
            push!(ph[i].strand2.window, 1)
        end

        # Second pass to find optimal break points and record info to phase
        @inbounds for w in 2:windows
            # get genotype vector spanning 2 windows
            start_prev = first(winranges[w - 1])
            start_curr = first(winranges[w])
            end_curr = last(winranges[w])
            Xwi = view(X, start_prev:end_curr, i)

            # previous and current haplotypes for both strands
            hap1_prev = haplotype1[i][w - 1]
            hap2_prev = haplotype2[i][w - 1]
            hap1_curr = haplotype1[i][w]
            hap2_curr = haplotype2[i][w]

            # find optimal breakpoint if there is one
            timers[id][16] += @elapsed _, bkpts, err = continue_haplotype(Xwi, 
                compressed_Hunique, w, (hap1_prev, hap2_prev),
                (hap1_curr, hap2_curr), phased=true, search_double_bkpts=true)

            timers[id][24] += @elapsed begin
                # record strand 1 info
                update_phase!(ph[i].strand1, compressed_Hunique, 
                    bkpts[1], hap1_prev, hap1_curr, Int32(w), 
                    start_prev, start_curr, end_curr, impute_untyped)
                # record strand 2 info
                update_phase!(ph[i].strand2, compressed_Hunique, 
                    bkpts[2], hap2_prev, hap2_curr, Int32(w), 
                    start_prev, start_curr, end_curr, impute_untyped)
            end
        end
        next!(pmeter) # update progress
    end
    return sum(timers) ./ Threads.nthreads()
end

function phase_fast_compressed!(
    ph::Vector{HaplotypeMosaicPair},
    X::AbstractMatrix{Union{Missing, T}},
    compressed_Hunique::CompressedHaplotypes,
    haplotype1::AbstractVector,
    haplotype2::AbstractVector,
    ) where T <: Real

    # declare some constants
    people = size(X, 2)
    snps = size(X, 1)
    haplotypes = nhaplotypes(compressed_Hunique)
    winranges = compressed_Hunique.X_window_range
    windows = length(haplotype1[1])

    # working arrays
    survivors1 = [Int32[] for _ in 1:Threads.nthreads()]
    survivors2 = [Int32[] for _ in 1:Threads.nthreads()]
    for id in 1:Threads.nthreads()
        sizehint!(survivors1[id], haplotypes)
        sizehint!(survivors2[id], haplotypes)
    end
    timers = [zeros(3*8) for _ in 1:Threads.nthreads()] # 8 for spacing
    pmeter = Progress(people, 5, "Phasing...")

    # phase person by person
    # for i in 1:people
    Threads.@threads for i in 1:people
        id = Threads.threadid()

        # First pass to phase each sample window-by-window
        timers[id][8] += @elapsed phase_sample!(haplotype1[i], haplotype2[i],
            compressed_Hunique, survivors1[id], survivors2[id])

        # record info for first window
        timers[id][24] += @elapsed begin
            h1 = haplotype1[i][1] # complete idx
            h2 = haplotype2[i][1] # complete idx
            push_Mosaic!(ph[i].strand1, (1, h1))
            push_Mosaic!(ph[i].strand2, (1, h2))
        end

        # Second pass to find optimal break points and record info to phase
        @inbounds for w in 2:windows
            # get genotype vector spanning 2 windows
            start_prev = first(winranges[w - 1])
            start_curr = first(winranges[w])
            end_curr = last(winranges[w])
            Xwi = view(X, start_prev:end_curr, i)

            # previous and current haplotypes for both strands
            hap1_prev = haplotype1[i][w - 1]
            hap2_prev = haplotype2[i][w - 1]
            hap1_curr = haplotype1[i][w]
            hap2_curr = haplotype2[i][w]

            # find optimal breakpoint if there is one
            timers[id][16] += @elapsed begin
                _, bkpts, err = continue_haplotype(Xwi, compressed_Hunique, w, 
                    (hap1_prev, hap2_prev), (hap1_curr, hap2_curr), 
                    phased=true, search_double_bkpts=true)
            end

            timers[id][24] += @elapsed begin
                # strand1 single stranded breakpoint
                if -1 < bkpts[1] < length(Xwi)
                    push_Mosaic!(ph[i].strand1, (start_prev + bkpts[1], 
                        hap1_curr))
                end
                # strand2 single stranded breakpoint
                if -1 < bkpts[2] < length(Xwi)
                    push_Mosaic!(ph[i].strand2, (start_prev + bkpts[2], 
                        hap2_curr))
                end
                # if not searching double bkpts, need to uncomment below
                # if bkpts[1] == -2
                #     push_Mosaic!(ph[i].strand1, (start_curr, hap1_curr))
                # end
                # if bkpts[2] == -2
                #     push_Mosaic!(ph[i].strand2, (start_curr, hap2_curr))
                # end
            end
        end
        next!(pmeter) # update progress
    end
    return sum(timers) ./ Threads.nthreads()
end
