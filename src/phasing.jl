"""
    phase(tgtfile, reffile; [outfile], [width], [flankwidth], [fast_method])

Phasing (haplotying) of `tgtfile` from a pool of haplotypes `reffile`
by sliding windows and saves result in `outfile`. Will also create an aligned 
reference `aligned.ref.vcf.gz` file matching `tgtfile` position by position. 

# Input
- `tgtfile`: VCF or PLINK file. VCF files should end in `.vcf` or `.vcf.gz`. PLINK files should exclude `.bim/.bed/.fam` suffixes but the trio must all be present in the directory.
- `reffile`: VCF file name. Should end in `.vcf` or `.vcf.gz`. 

# Optional Inputs
- `outfile`: output filename ending in `.vcf.gz` or `.vcf`. Output genotypes will have no missing data.
- `impute`: If `true`, untyped SNPs will be imputed, otherwise only missing snps in `tgtfile` will be imputed.  (default `false`)
- `width`: number of SNPs (markers) in each haplotype window. (default `2000`)
- `flankwidth`: Number of SNPs flanking the sliding window (defaults to 10% of `width`)
- `typed_snps_threshold`: Number of typed SNPs required in each window. Below this threshold means optimal haplotype pair will be copied from the nearest window with enough typed SNPs
"""
function phase(
    tgtfile::AbstractString,
    reffile::AbstractString;
    outfile::AbstractString = "imputed." * tgtfile,
    impute::Bool = true, #TODO for impute=false
    width::Int = 2000,
    min_typed_snps = 100, 
    fast_method::Bool = false,
    unique_only::Bool = false
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
        # filter for unique haplotypes if not already done 
        compressed_Hunique = compress_haplotypes(reffile, "compressed." * reffile, width, dims=2, flankwidth = 0)
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
    ref_snps = compressed_Hunique.snps
    windows = floor(Int, ref_snps / width)

    #
    # compute redundant haplotype sets
    #
    calculate_happairs_start = time()
    pmeter = Progress(windows, 5, "Computing optimal haplotype pairs...")
    redundant_haplotypes = [[Tuple{Int, Int}[] for i in 1:windows] for j in 1:people]
    typed_snps = Vector{Vector{Int}}(undef, windows) #tracks index for typed snps in each window
    mutex = Threads.SpinLock()
    Threads.@threads for w in 1:windows
        # match target and ref file by snp position
        Threads.lock(mutex)
        cur_range = compressed_Hunique.CWrange[w]
        Hw_pos = compressed_Hunique.pos[cur_range]
        XtoH_idx = indexin(X_pos, Hw_pos) # X_pos[i] == Hw_pos[XtoH_idx[i]]
        XtoH_rm_nothing = Base.filter(!isnothing, XtoH_idx) # delete snps not in ref panel
        Xw_aligned = X[findall(!isnothing, XtoH_idx), :]
        Hw_aligned = compressed_Hunique[w].uniqueH[XtoH_rm_nothing, :]
        typed_snps[w] = XtoH_rm_nothing # save typed snps index for current window
        Threads.unlock(mutex)

        # Skip windows with too few typed SNPs
        if length(XtoH_rm_nothing) < min_typed_snps
            for k in 1:people
                push!(redundant_haplotypes[k][w], (-1, -1))
            end
            continue
        end

        # computational routine (TODO: preallocate all internal matrices here)
        happairs, hapscore = haplopair(Xw_aligned, Hw_aligned)
        
        # convert happairs (which index off unique haplotypes) to indices of full haplotype pool, and find all matching happairs
        compute_redundant_haplotypes!(redundant_haplotypes, compressed_Hunique, happairs, w)

        # update progress
        next!(pmeter)
    end
    calculate_happairs_time = time() - calculate_happairs_start

    #
    # phasing (haplotyping) step
    #
    # offset = (chunks - 1) * snps_per_chunk
    phase_start = time()
    ph = [HaplotypeMosaicPair(ref_snps) for i in 1:people]
    phase!(ph, X, compressed_Hunique, typed_snps, redundant_haplotypes, X_pos) # phase by dynamic programming + breakpoint search
    phase_time = time() - phase_start

    #
    # impute step
    #
    impute_start = time()
    H_pos = compressed_Hunique.pos
    XtoH_idx = indexin(X_pos, H_pos) # X_pos[i] == H_pos[XtoH_idx[i]]
    XtoH_rm_nothing = Base.filter(!isnothing, XtoH_idx)
    X_aligned = any(isnothing.(XtoH_idx)) ? X[findall(!isnothing, XtoH_idx), :] : X
    X_full = Matrix{Union{Missing, UInt8}}(missing, ref_snps, people)
    copyto!(@view(X_full[XtoH_rm_nothing, :]), X_aligned)
    impute!(X_full, compressed_Hunique, ph, outfile, X_sampleID) # imputes X_full and writes to file
    impute_time = time() - impute_start

    println("Data import time                    = ", round(import_data_time, sigdigits=6), " seconds")
    println("Computing haplotype pair time       = ", round(calculate_happairs_time, sigdigits=6), " seconds")
    println("Phasing by dynamic programming time = ", round(phase_time, sigdigits=6), " seconds")
    println("Imputing time                       = ", round(impute_time, sigdigits=6), " seconds")

    return ph, redundant_haplotypes
end

"""
    phase!(X, H, width=400, verbose=true)

Phasing (haplotying) of genotype matrix `X` from a pool of haplotypes `H`
by dynamic programming. A precomputed, window-by-window haplotype pairs is assumed. 

# Input
* `ph`: A vector of `HaplotypeMosaicPair` keeping track of each person's phase information.
* `X`: `p x n` matrix with missing values. Each column is genotypes of an individual.
* `compressed_Hunique`: A `CompressedHaplotypes` keeping track of unique haplotypes for each window and some other information
* `typed_snps`: `typed_snps[w]` are indices of typed SNPs in window `w`
* `hapset`: Vector of optimal haplotype pairs across windows. The haplotype pairs are indices to the full haplotype set and NOT the compressed haplotypes
* `chunk_offset`: Shifts SNPs if a chromosome had been chunked. (not currently implemented)
"""
function phase!(
    ph::Vector{HaplotypeMosaicPair},
    X::AbstractMatrix{Union{Missing, T}},
    compressed_Hunique::CompressedHaplotypes,
    typed_snps::Vector{Vector{Int}},
    hapset::Vector{Vector{Vector{Tuple{Int, Int}}}},
    X_pos::Vector{Int};
    chunk_offset::Int = 0,
    ) where T <: Real

    # declare some constants
    people = size(X, 2)
    haplotypes = nhaplotypes(compressed_Hunique)
    snps = compressed_Hunique.snps
    width = compressed_Hunique.width
    windows = floor(Int, snps / width)
    H_pos = compressed_Hunique.pos
    last_window_width = snps - (windows - 1) * width 

    # allocate working arrays
    sol_path = [Vector{Tuple{Int, Int}}(undef, windows) for i in 1:Threads.nthreads()]
    nxt_pair = [[Int[] for i in 1:windows] for i in 1:Threads.nthreads()]
    tree_err = [[Float64[] for i in 1:windows] for i in 1:Threads.nthreads()]
    HtoX_idx = indexin(H_pos, X_pos)
    XtoH_idx = indexin(X_pos, H_pos)
    pmeter   = Progress(people, 5, "Merging breakpoints...")

    # loop over each person
    # first  1/3: ((w - 2) * width + 1):((w - 1) * width)
    # middle 1/3: ((w - 1) * width + 1):(      w * width)
    # last   1/3: (      w * width + 1):((w + 1) * width)
    for i in 1:people
        # first find optimal haplotype pair in each window using dynamic programming
        id = Threads.threadid()
        connect_happairs!(sol_path[id], nxt_pair[id], tree_err[id], hapset[i], λ = 1.0)

        # phase first window 
        k, l = sol_path[id][1][1], sol_path[id][1][2] # complete haplotype index
        h1 = complete_idx_to_unique_idx(k, 1, compressed_Hunique)
        h2 = complete_idx_to_unique_idx(l, 1, compressed_Hunique)
        push!(ph[i].strand1.start, 1 + chunk_offset)
        push!(ph[i].strand1.window, 1) 
        push!(ph[i].strand1.haplotypelabel, h1)
        push!(ph[i].strand2.start, 1 + chunk_offset)
        push!(ph[i].strand2.window, 1)
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
        #     h1 = complete_idx_to_unique_idx(k, w, compressed_Hunique)
        #     h2 = complete_idx_to_unique_idx(l, w, compressed_Hunique)

        #     push!(ph[i].strand1.start, chunk_offset + (w - 1) * width + 1)
        #     push!(ph[i].strand1.haplotypelabel, h1)
        #     push!(ph[i].strand1.window, w)
        #     push!(ph[i].strand2.start, chunk_offset + (w - 1) * width + 1)
        #     push!(ph[i].strand2.haplotypelabel, h2)
        #     push!(ph[i].strand2.window, w)
        # end

        # search breakpoints then record result into ph
        for w in 2:windows
            # get imputation target range = 2 windows of H, tracking untyped snps
            Hw_start  = (w - 2) * width + 1
            Hw_mid    = (w - 1) * width + 1
            Hw_end    = (w == windows ? snps : w * width)
            Xwi_start = HtoX_idx[something(findnext(!isnothing, HtoX_idx, Hw_start))]
            Xwi_end   = HtoX_idx[something(findprev(!isnothing, HtoX_idx, Hw_end))]
            Xwi = view(X, Xwi_start:Xwi_end, i)

            # find optimal breakpoint
            sol_path[id][w], bkpts = continue_haplotype(Xwi, compressed_Hunique, 
                typed_snps, w, sol_path[id][w - 1], sol_path[id][w])

            # record strand 1 info
            update_phase!(ph[i].strand1, compressed_Hunique, bkpts[1], sol_path[id][w - 1][1], 
                sol_path[id][w][1], w, width, chunk_offset, Hw_start, Hw_mid, Hw_end, 
                HtoX_idx, XtoH_idx, Xwi_start)
            # record strand 2 info
            update_phase!(ph[i].strand2, compressed_Hunique, bkpts[2], sol_path[id][w - 1][2], 
                sol_path[id][w][2], w, width, chunk_offset, Hw_start, Hw_mid, Hw_end, 
                HtoX_idx, XtoH_idx, Xwi_start)
        end

        # update progress
        next!(pmeter)
    end
end

"""
Helper function for updating phase information after breakpoints have been identified
between windows `w - 1` and `w`. Every window have 0 or 1 breakpoint.
"""
function update_phase!(ph::HaplotypeMosaic, compressed_Hunique::CompressedHaplotypes,
    bkpt::Int, hap_prev::Int, hap_curr::Int, w::Int, width::Int, chunk_offset::Int,
    Hw_start::Int, Hw_mid::Int, Hw_end::Int, HtoX_idx::AbstractVector, 
    XtoH_idx::AbstractVector, Xwi_start::Int)

    # no breakpoints
    if bkpt == -1
        h = complete_idx_to_unique_idx(hap_curr, w, compressed_Hunique)
        push!(ph.start, chunk_offset + (w - 1) * width + 1)
        push!(ph.haplotypelabel, h)
        push!(ph.window, w)
        return nothing
    end

    # convert bkpt (in terms of X index) to index in H
    X_bkpt_end = Xwi_start + bkpt
    H_bkpt_pos = XtoH_idx[X_bkpt_end]

    if H_bkpt_pos >= Hw_end
        # previous window dominates current window 
        h = complete_idx_to_unique_idx(hap_prev, w, compressed_Hunique)
        push!(ph.start, chunk_offset + (w - 1) * width + 1)
        push!(ph.haplotypelabel, h)
        push!(ph.window, w)
    elseif Hw_mid <= H_bkpt_pos < Hw_end
        # previous window extends to current window 
        h1 = complete_idx_to_unique_idx(hap_prev, w, compressed_Hunique)
        push!(ph.start, chunk_offset + (w - 1) * width + 1)
        push!(ph.haplotypelabel, h1)
        push!(ph.window, w)
        # 2nd part of current window
        h2 = complete_idx_to_unique_idx(hap_curr, w, compressed_Hunique)
        X_bkpt_end = Xwi_start + bkpt
        H_bkpt_pos = XtoH_idx[X_bkpt_end]
        push!(ph.start, chunk_offset + H_bkpt_pos)
        push!(ph.haplotypelabel, h2)
        push!(ph.window, w)
    elseif H_bkpt_pos < Hw_mid
        # current window extends to previous window
        h1 = complete_idx_to_unique_idx(hap_curr, w - 1, compressed_Hunique)
        X_bkpt_end = Xwi_start + bkpt
        H_bkpt_pos = XtoH_idx[X_bkpt_end]
        push!(ph.start, chunk_offset + H_bkpt_pos)
        push!(ph.haplotypelabel, h1)
        push!(ph.window, w - 1)
        # update current window
        h2 = complete_idx_to_unique_idx(hap_curr, w, compressed_Hunique)
        push!(ph.start, chunk_offset + (w - 1) * width + 1)
        push!(ph.haplotypelabel, h2)
        push!(ph.window, w)
    else
        # println("H_bkpt_pos = $H_bkpt_pos, Hw_start=$Hw_start, Hw_mid=$Hw_mid, Hw_end=$Hw_end ")
        error("update_phase!: bkpt does not satisfy -1 <= bkpt <= 2width! Shouldn't be possible")
    end

    return nothing
end

function phase_fast!(
    ph::Vector{HaplotypeMosaicPair},
    X::AbstractMatrix{Union{Missing, T}},
    H::AbstractMatrix,
    hapset::Vector{OptimalHaplotypeSet};
    width::Int = 400,
    flankwidth::Int = round(Int, 0.1width),
    chunk_offset::Int = 0,
    Xtrue::Union{AbstractMatrix, Nothing} = nothing, # for testing
    ) where T <: Real

    # declare some constants
    snps, people = size(X)
    haplotypes = size(H, 2)
    windows = floor(Int, snps / width)

    # allocate working arrays
    haplo_chain = ([copy(hapset[i].strand1[1]) for i in 1:people], [copy(hapset[i].strand2[1]) for i in 1:people])
    chain_next  = (BitVector(undef, haplotypes), BitVector(undef, haplotypes))
    window_span = (ones(Int, people), ones(Int, people))
    pmeter      = Progress(people, 1, "Intersecting haplotypes...")

    # begin intersecting haplotypes window by window
    @inbounds for i in 1:people
        for w in 2:windows
            # Decide whether to cross over based on the larger intersection
            # A   B      A   B
            # |   |  or    X
            # C   D      C   D
            chain_next[1] .= haplo_chain[1][i] .& hapset[i].strand1[w] # not crossing over
            chain_next[2] .= haplo_chain[1][i] .& hapset[i].strand2[w] # crossing over
            AC = sum(chain_next[1])
            AD = sum(chain_next[2])
            chain_next[1] .= haplo_chain[2][i] .& hapset[i].strand1[w] # crossing over
            chain_next[2] .= haplo_chain[2][i] .& hapset[i].strand2[w] # not crossing over
            BC = sum(chain_next[1])
            BD = sum(chain_next[2])
            if AC + BD < AD + BC
                hapset[i].strand1[w], hapset[i].strand2[w] = hapset[i].strand2[w], hapset[i].strand1[w]
            end

            # intersect all surviving haplotypes with next window
            chain_next[1] .= haplo_chain[1][i] .& hapset[i].strand1[w]
            chain_next[2] .= haplo_chain[2][i] .& hapset[i].strand2[w]

            # strand 1 becomes empty
            if sum(chain_next[1]) == 0
                # delete all nonmatching haplotypes in previous windows
                for ww in (w - window_span[1][i]):(w - 1)
                    hapset[i].strand1[ww] .= haplo_chain[1][i]
                end

                # reset counters and storage
                haplo_chain[1][i] .= hapset[i].strand1[w]
                window_span[1][i] = 1
            else
                haplo_chain[1][i] .= chain_next[1]
                window_span[1][i] += 1
            end

            # strand 2 becomes empty
            if sum(chain_next[2]) == 0
                # delete all nonmatching haplotypes in previous windows
                for ww in (w - window_span[2][i]):(w - 1)
                    hapset[i].strand2[ww] .= haplo_chain[2][i]
                end

                # reset counters and storage
                haplo_chain[2][i] .= hapset[i].strand2[w]
                window_span[2][i] = 1
            else
                haplo_chain[2][i] .= chain_next[2]
                window_span[2][i] += 1
            end
        end
        next!(pmeter) #update progress
    end

    # get rid of redundant haplotypes in last few windows separately, since intersection may not become empty
    for i in 1:people
        for ww in (windows - window_span[1][i] + 1):windows
            hapset[i].strand1[ww] .= haplo_chain[1][i]
        end

        for ww in (windows - window_span[2][i] + 1):windows
            hapset[i].strand2[ww] .= haplo_chain[2][i]
        end
    end

    # phase window 1
    for i in 1:people
        hap1 = findfirst(hapset[i].strand1[1]) :: Int64
        hap2 = findfirst(hapset[i].strand2[1]) :: Int64
        push!(ph[i].strand1.start, 1 + chunk_offset)
        push!(ph[i].strand1.haplotypelabel, hap1)
        push!(ph[i].strand2.start, 1 + chunk_offset)
        push!(ph[i].strand2.haplotypelabel, hap2)
    end

    # find optimal break points and record info to phase
    # first  1/3: ((w - 2) * width + 1):((w - 1) * width)
    # middle 1/3: ((w - 1) * width + 1):(      w * width)
    # last   1/3: (      w * width + 1):((w + 1) * width)
    pmeter = Progress(people, 5, "Merging breakpoints...")
    strand1_intersect = [copy(chain_next[1]) for i in 1:Threads.nthreads()]
    strand2_intersect = [copy(chain_next[2]) for i in 1:Threads.nthreads()]
    Threads.@threads for i in 1:people
        id = Threads.threadid()
        for w in 2:windows
            w_start = max(ph[i].strand1.start[end], ph[i].strand2.start[end], ((w - 2) * width + 1))
            Hi = view(H, w_start:(w * width), :)
            Xi = view(X, w_start:(w * width), i)
            strand1_intersect[id] .= hapset[i].strand1[w - 1] .& hapset[i].strand1[w]
            strand2_intersect[id] .= hapset[i].strand2[w - 1] .& hapset[i].strand2[w]
            # double haplotype switch
            if sum(strand1_intersect[id]) == sum(strand2_intersect[id]) == 0
                s1_prev = ph[i].strand1.haplotypelabel[end]
                s2_prev = ph[i].strand2.haplotypelabel[end]
                # search breakpoints when choosing first pair
                s1_next = findfirst(hapset[i].strand1[w]) :: Int64
                s2_next = findfirst(hapset[i].strand2[w]) :: Int64
                bkpt, err_optim = search_breakpoint(Xi, Hi, (s1_prev, s1_next), (s2_prev, s2_next))
                # record info into phase
                push!(ph[i].strand1.start, chunk_offset + w_start + bkpt[1])
                push!(ph[i].strand2.start, chunk_offset + w_start + bkpt[2])
                push!(ph[i].strand1.haplotypelabel, s1_next)
                push!(ph[i].strand2.haplotypelabel, s2_next)
            else
                # single haplotype switch
                if sum(strand1_intersect[id]) == 0
                    # search breakpoints among all possible haplotypes
                    s1_prev = ph[i].strand1.haplotypelabel[end]
                    s1_win_next = findall(hapset[i].strand1[w])
                    s2_win_next = findall(hapset[i].strand2[w])
                    best_bktp = 0
                    best_err  = typemax(Int)
                    best_s1_next = 0
                    for s1_next in s1_win_next, s2_next in s2_win_next
                        bkpt, err_optim = search_breakpoint(Xi, Hi, s2_next, (s1_prev, s1_next))
                        if err_optim < best_err
                            best_bktp, best_err, best_s1_next = bkpt, err_optim, s1_next
                        end
                    end
                    # record info into phase
                    push!(ph[i].strand1.start, chunk_offset + w_start + best_bktp)
                    push!(ph[i].strand1.haplotypelabel, best_s1_next)
                end

                # single haplotype switch
                if sum(strand2_intersect[id]) == 0
                    # search breakpoints among all possible haplotypes
                    s2_prev = ph[i].strand2.haplotypelabel[end]
                    s2_win_next = findall(hapset[i].strand2[w])
                    s1_win_next = findall(hapset[i].strand1[w])
                    best_bktp = 0
                    best_err  = typemax(Int)
                    best_s2_next = 0
                    for s2_next in s2_win_next, s1_next in s1_win_next
                        bkpt, err_optim = search_breakpoint(Xi, Hi, s1_next, (s2_prev, s2_next))
                        if err_optim < best_err
                            best_bktp, best_err, best_s2_next = bkpt, err_optim, s2_next
                        end
                    end
                    # record info into phase
                    push!(ph[i].strand2.start, chunk_offset + w_start + best_bktp)
                    push!(ph[i].strand2.haplotypelabel, best_s2_next)
                end
            end
        end
        next!(pmeter) #update progress
    end
end


"""
    phase_dp_fast!(X, H, width=400, verbose=true)
Phasing (haplotying) of genotype matrix `X` from a pool of haplotypes `H`
by dynamic programming. A precomputed, window-by-window haplotype pairs is assumed. 
# Input
* `X`: `p x n` matrix with missing values. Each column is genotypes of an individual.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `width`: width of the sliding window.
* `verbose`: display algorithmic information.
"""
function phase_unique_only!(
    ph::Vector{HaplotypeMosaicPair},
    X::AbstractMatrix{Union{Missing, T}},
    H::AbstractMatrix,
    hapset::Vector{Vector{Vector{Tuple{Int, Int}}}};
    width::Int = 400,
    flankwidth::Int = round(Int, 0.1width),
    chunk_offset::Int = 0,
    Xtrue::Union{AbstractMatrix, Nothing} = nothing, # for testing
    ) where T <: Real

    # declare some constants
    snps, people = size(X)
    haplotypes = size(H, 2)
    windows = floor(Int, snps / width)
    pmeter = Progress(people, 5, "Imputing typed markers...")

    # loop over each person
    # first  1/3: ((w - 2) * width + 1):((w - 1) * width)
    # middle 1/3: ((w - 1) * width + 1):(      w * width)
    # last   1/3: (      w * width + 1):((w + 1) * width)
    Threads.@threads for i in 1:people
        # phase first window 
        push!(ph[i].strand1.start, 1 + chunk_offset)
        push!(ph[i].strand1.haplotypelabel, hapset[i][1][1][1]) # i'th person, 1st window, 1st haplotype pair, 1st haplotype
        push!(ph[i].strand2.start, 1 + chunk_offset)
        push!(ph[i].strand2.haplotypelabel, hapset[i][1][1][2])

        # phase without searching for breakpoints
        for w in 2:windows
            u, j = hapset[i][w - 1][1] # 1st haplotype pair in previous window for person i
            k, l = hapset[i][w][1]     # 1st haplotype pair in current window for person i

            # switch current window's pair order if 1 or 2 haplotype match
            if (u == l && j == k) || (j == k && u ≠ l) || (u == l && j ≠ k)
                k, l = l, k 
            end

            push!(ph[i].strand1.start, chunk_offset + (w - 1) * width + 1)
            push!(ph[i].strand1.haplotypelabel, k)
            push!(ph[i].strand2.start, chunk_offset + (w - 1) * width + 1)
            push!(ph[i].strand2.haplotypelabel, l)
        end

        # phase middle windows
        # for w in 2:(windows - 1)
        #     Xwi = view(X, ((w - 2) * width + 1):(w * width), i)
        #     Hw  = view(H, ((w - 2) * width + 1):(w * width), :)
        #     hapset[i][w], bkpts = continue_haplotype(Xwi, Hw, hapset[i][w - 1], hapset[i][w])

        #     # strand 1
        #     if bkpts[1] > -1 && bkpts[1] < 2width
        #         push!(ph[i].strand1.start, chunk_offset + (w - 2) * width + 1 + bkpts[1])
        #         push!(ph[i].strand1.haplotypelabel, hapset[i][w][1])
        #     end
        #     # strand 2
        #     if bkpts[2] > -1 && bkpts[2] < 2width
        #         push!(ph[i].strand2.start, chunk_offset + (w - 2) * width + 1 + bkpts[2])
        #         push!(ph[i].strand2.haplotypelabel, hapset[i][w][2])
        #     end
        # end

        # # phase last window
        # Xwi = view(X, ((windows - 2) * width + 1):snps, i)
        # Hw  = view(H, ((windows - 2) * width + 1):snps, :)
        # hapset[i][windows], bkpts = continue_haplotype(Xwi, Hw, hapset[i][windows - 1], hapset[i][windows])
        # # strand 1
        # if bkpts[1] > -1 && bkpts[1] < 2width
        #     push!(ph[i].strand1.start, chunk_offset + (windows - 2) * width + 1 + bkpts[1])
        #     push!(ph[i].strand1.haplotypelabel, hapset[i][windows][1])
        # end
        # # strand 2
        # if bkpts[2] > -1 && bkpts[2] < 2width
        #     push!(ph[i].strand2.start, chunk_offset + (windows - 2) * width + 1 + bkpts[2])
        #     push!(ph[i].strand2.haplotypelabel, hapset[i][windows][2])
        # end

        # update progress
        next!(pmeter) 
    end
end
