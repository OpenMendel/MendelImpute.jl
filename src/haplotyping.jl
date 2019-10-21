"""
    phase(X, H, width=400, verbose=true)

Phasing (haplotying) of genotype matrix `X` from a pool of haplotypes `H`
by sliding windows.

# Input
* `X`: `p x n` matrix with missing values. Each column is genotypes of an individual.
* `H`: `p x d` haplotype matrix. Each column is a haplotype.
* `width`: width of the sliding window.
* `verbose`: display algorithmic information.
"""
function phase(
    X::AbstractMatrix{Union{Missing, T}},
    H::AbstractMatrix{T},
    width::Int    = 400,
    verbose::Bool = true
    ) where T <: Real

    people, snps, haplotypes = size(X, 2), size(X, 1), size(H, 2)
    # allocate working arrays
    M        = zeros(T, haplotypes, haplotypes)
    N        = zeros(T,     people, haplotypes)
    happair  = ones(Int, people), ones(Int, people)
    hapscore = zeros(T, people)
    phase    = [HaplotypeMosaicPair(snps) for i in 1:people]

    # no need for sliding window
    if snps ≤ 3width
        haploimpute!(X, H, M, N, happair, hapscore)
        for i in 1:people
            push!(phase[i].strand1.start, 1)
            push!(phase[i].strand1.haplotypelabel, happair[1][i])
            push!(phase[i].strand2.start, 1)
            push!(phase[i].strand2.haplotypelabel, happair[2][i])
        end
        return phase
    end

    # allocate working arrays
    Xwork = X[1:3width, :]
    Xwork_float = zeros(T, size(Xwork))
    Hwork = view(H, 1:3width, :)
    # Hwork, (pp, dd) = unique_haplotypes(H, ((windows - 2) * width + 1):snps)
    happair_prev = deepcopy(happair)

    # number of windows
    windows = floor(Int, snps / width)

    # phase and impute window 1
    verbose && println("Imputing SNPs 1:$width")
    haploimpute!(Xwork, Hwork, M, N, happair, hapscore, Xfloat=Xwork_float)
    for i in 1:people
        push!(phase[i].strand1.start, 1)
        push!(phase[i].strand1.haplotypelabel, happair[1][i])
        push!(phase[i].strand2.start, 1)
        push!(phase[i].strand2.haplotypelabel, happair[2][i])
    end

    # first  1/3: ((w - 2) * width + 1):((w - 1) * width)
    # middle 1/3: ((w - 1) * width + 1):(      w * width)
    # last   1/3: (      w * width + 1):((w + 1) * width)
    for w in 2:(windows - 1)
        if verbose
            println("Imputing SNPs $((w - 1) * width + 1):$(w * width)")
        end

        # sync Xwork and Hwork with original data
        Hwork = view(H, ((w - 2) * width + 1):((w + 1) * width), :)
        # Hwork, (pp, dd) = unique_haplotypes(H, ((w - 2) * width + 1):((w + 1) * width))
        copyto!(Xwork, view(X, ((w - 2) * width + 1):((w + 1) * width), :))

        # phase current window
        copyto!(happair_prev[1], happair[1])
        copyto!(happair_prev[2], happair[2])
        haploimpute!(Xwork, Hwork, M, N, happair, hapscore, Xfloat=Xwork_float)
        # haploimpute!(Xwork, Hwork, view(M, 1:dd, 1:dd), view(N, :, 1:dd), happair, hapscore, Xfloat=Xwork_float)

        # find optimal break points and record info into phase
        Hw12 = view(Hwork, 1:2width, :)
        for i in 1:people
            Xi = view(Xwork, 1:2width, i)
            (happair[1][i], happair[2][i]), bkpts =
                continue_haplotype(Xi, Hw12,
                (happair_prev[1][i], happair_prev[2][i]),
                (     happair[1][i],      happair[2][i]))
            # strand 1
            if bkpts[1] > -1 && bkpts[1] < 2width
                push!(phase[i].strand1.start, (w - 2) * width + 1 + bkpts[1])
                push!(phase[i].strand1.haplotypelabel, happair[1][i])
            end
            # strand 2
            if bkpts[2] > -1 && bkpts[2] < 2width
                push!(phase[i].strand2.start, (w - 2) * width + 1 + bkpts[2])
                push!(phase[i].strand2.haplotypelabel, happair[2][i])
            end
            # # for debug
            if verbose == true && i == 1
                println("happair = ($(happair[1][i]), $(happair[2][i]))")
                println("bkpts = $bkpts")
            end
        end
    end

    # Hua's code without searching breakpoints
    # for w in 2:(windows-1)
    #     if verbose
    #         println("Imputing SNPs $((w - 1) * width + 1):$(w * width)")
    #     end

    #     # sync Xwork and Hwork with original data
    #     Hwork = view(H, ((w - 2) * width + 1):((w + 1) * width), :)
    #     # Hwork, (pp, dd) = unique_haplotypes(H, ((w - 2) * width + 1):((w + 1) * width))
    #     copyto!(Xwork, view(X, ((w - 2) * width + 1):((w + 1) * width), :))

    #     # phase current window
    #     copyto!(happair_prev[1], happair[1])
    #     copyto!(happair_prev[2], happair[2])
    #     haploimpute!(Xwork, Hwork, M, N, happair, hapscore, Xfloat=Xwork_float)
    #     # haploimpute!(Xwork, Hwork, view(M, 1:dd, 1:dd), view(N, :, 1:dd), happair, hapscore, Xfloat=Xwork_float)

    #     # record info into phase
    #     for i in 1:people
    #         push!(phase[i].strand1.start, (w - 2) * width + 1)
    #         push!(phase[i].strand1.haplotypelabel, happair[1][i])
    #         push!(phase[i].strand2.start, (w - 2) * width + 1)
    #         push!(phase[i].strand2.haplotypelabel, happair[2][i])
    #     end
    # end

    # phase last window
    if verbose
        println("Imputing SNPs $((windows - 1) * width + 1):$snps")
    end
    Xwork = X[((windows - 2) * width + 1):snps, :]
    Hwork = view(H, ((windows - 2) * width + 1):snps, :)
    # Hwork, (pp, dd) = unique_haplotypes(H, ((windows - 2) * width + 1):snps)
    copyto!(happair_prev[1], happair[1])
    copyto!(happair_prev[2], happair[2])
    # haploimpute!(Xwork, Hwork, view(M, 1:dd, 1:dd), view(N, :, dd), happair, hapscore)
    haploimpute!(Xwork, Hwork, M, N, happair, hapscore)

    # # find optimal break points and record info to phase
    for i in 1:people
        (happair[1][i], happair[2][i]), bkpts =
        continue_haplotype(Xwork[:, i], Hwork,
            (happair_prev[1][i], happair_prev[2][i]),
            (happair[1][i], happair[2][i]))
        # strand 1
        if bkpts[1] > -1 && bkpts[1] < 2width
            push!(phase[i].strand1.start, (windows - 2) * width + 1 + bkpts[1])
            push!(phase[i].strand1.haplotypelabel, happair[1][i])
        end
        # strand 2
        if bkpts[2] > -1 && bkpts[2] < 2width
            push!(phase[i].strand2.start, (windows - 2) * width + 1 + bkpts[2])
            push!(phase[i].strand2.haplotypelabel, happair[2][i])
        end
    end

    return phase
end

function phase2(
    X::AbstractMatrix{Union{Missing, T}},
    H::AbstractMatrix{T};
    width::Int    = 128,
    verbose::Bool = true
    ) where T <: Real

    # problem dimensions
    snps, people = size(X)

    # number of windows
    windows = floor(Int, snps / width)

    # get redundant haplotype sets. 
    hapset = compute_optimal_halotype_set(X, H, width=width, verbose=verbose)

    # allocate working arrays
    phase = [HaplotypeMosaicPair(snps) for i in 1:people]
    store = ([copy(hapset.strand1[1, i]) for i in 1:people], [copy(hapset.strand2[1, i]) for i in 1:people])
    window_span = (ones(Int, people), ones(Int, people))

    # TODO: parallel computing
    # TODO: replace `intersect` and `intersect!` with fast set intersection using bisection/seesaw search
    @inbounds for i in 1:people, w in 2:windows

        # decide how to concatenate next 2 windows to previous windows based on the larger intersection
        A = intersect(store[1][i], hapset.strand1[w, i])
        B = intersect(store[1][i], hapset.strand2[w, i])
        if length(A) >= length(B)
            # no need to cross over
            a = A 
        else
            # cross over
            a = B
            hapset.strand1.p[w, i], hapset.strand2.p[w, i] = hapset.strand2.p[w, i], hapset.strand1.p[w, i]
        end
        b = intersect(store[2][i], hapset.strand2[w, i])

        # strand 1
        if isempty(a)
            # delete all nonmatching haplotypes in previous windows
            for ww in (w - window_span[1][i]):(w - 1)
                hapset.strand1.p[ww, i] = copy(store[1][i]) 
            end

            # update counters and storage
            store[1][i] = copy(hapset.strand1[w, i])
            window_span[1][i] = 1
        else
            intersect!(store[1][i], hapset.strand1[w, i])
            window_span[1][i] += 1
        end

        # strand 2
        if isempty(b)
            # delete all nonmatching haplotypes in previous windows
            for ww in (w - window_span[2][i]):(w - 1)
                hapset.strand2.p[ww, i] = copy(store[2][i]) 
            end

            # update counters and storage
            store[2][i] = copy(hapset.strand2[w, i])
            window_span[2][i] = 1
        else
            intersect!(store[2][i], hapset.strand2[w, i])
            window_span[2][i] += 1
        end
    end

    # TODO: there's a bug in computing redundant haplotypes since last window never agrees with 2nd to last window
    # handle last few windows separately, since they may not hit the isempty command
    for i in 1:people
        for ww in (windows - window_span[1][i] + 1):windows
            hapset.strand1.p[ww, i] = copy(store[1][i]) 
        end

        for ww in (windows - window_span[2][i] + 1):windows
            hapset.strand2.p[ww, i] = copy(store[2][i]) 
        end
    end

    # phase window 1
    for i in 1:people
        push!(phase[i].strand1.start, 1)
        push!(phase[i].strand1.haplotypelabel, first(hapset.strand1[1, i]))
        push!(phase[i].strand2.start, 1)
        push!(phase[i].strand2.haplotypelabel, first(hapset.strand2[1, i]))
    end

    #phase window by window without checking breakpoints
    for i in 1:people, w in 2:windows
        hap1 = first(hapset.strand1[w, i])
        hap2 = first(hapset.strand2[w, i])

        # strand 1
        push!(phase[i].strand1.start, (w - 1) * width + 1)
        push!(phase[i].strand1.haplotypelabel, hap1)

        # strand 2
        push!(phase[i].strand2.start, (w - 1) * width + 1)
        push!(phase[i].strand2.haplotypelabel, hap2)
    end

    # find optimal break points and record info to phase. 
    # store = ([copy(hapset.strand1[1, i]) for i in 1:people], [copy(hapset.strand2[1, i]) for i in 1:people])
    # for i in 1:people, w in 2:windows
        
    #     a = intersect(store[1][i], hapset.strand1[w, i])
    #     b = intersect(store[2][i], hapset.strand2[w, i])

    #     if isempty(a)
    #         # search breakpoints
    #         Xi = view(X, ((w - 2) * width + 1):(w * width), i)
    #         Hi = view(H, ((w - 2) * width + 1):(w * width), :)
    #         prev_and_cur_haplotypes = (hapset.strand1[w - 1, i], hapset.strand1[w, i])
    #         bkpt, hap, err_optim = search_breakpoint(Xi, Hi, hapset.strand2[w, i], prev_and_cur_haplotypes)

    #         # record info into phase
    #         push!(phase[i].strand1.start, (w - 2) * width + 1 + bkpt)
    #         push!(phase[i].strand1.haplotypelabel, hap)

    #         # update storage
    #         store[1][i] = copy(hapset.strand1[w, i])
    #     end

    #     if isempty(b)
    #         # search breakpoints
    #         Xi = view(X, ((w - 2) * width + 1):(w * width), i)
    #         Hi = view(H, ((w - 2) * width + 1):(w * width), :)
    #         prev_and_cur_haplotypes = (hapset.strand2[w - 1, i], hapset.strand2[w, i])
    #         bkpt, hap, err_optim = search_breakpoint(Xi, Hi, hapset.strand1[w, i], prev_and_cur_haplotypes)

    #         # record info into phase
    #         push!(phase[i].strand2.start, (w - 2) * width + 1 + bkpt)
    #         push!(phase[i].strand2.haplotypelabel, hap)

    #         # update storage
    #         store[2][i] = copy(hapset.strand2[w, i])
    #     end
    # end

    # finally, fill in missing entries of X
    # impute!(X, H, phase)
    impute2!(X, H, phase)

    return hapset
    # return phase, hapset, bkpts
end

function impute!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    phase::Vector{HaplotypeMosaicPair}
    )

    fill!(X, 0)
    # loop over individuals
    for i in 1:size(X, 2)
        for s in 1:(length(phase[i].strand1.start) - 1)
            idx = phase[i].strand1.start[s]:(phase[i].strand1.start[s + 1] - 1)
            X[idx, i] = H[idx, phase[i].strand1.haplotypelabel[s]]
        end
        idx = phase[i].strand1.start[end]:phase[i].strand1.length
        X[idx, i] = H[idx, phase[i].strand1.haplotypelabel[end]]
        for s in 1:(length(phase[i].strand2.start) - 1)
            idx = phase[i].strand2.start[s]:(phase[i].strand2.start[s + 1] - 1)
            X[idx, i] += H[idx, phase[i].strand2.haplotypelabel[s]]
        end
        idx = phase[i].strand2.start[end]:phase[i].strand2.length
        X[idx, i] += H[idx, phase[i].strand2.haplotypelabel[end]]
    end
end

function impute2!(
    X::AbstractMatrix,
    H::AbstractMatrix,
    phase::Vector{HaplotypeMosaicPair}
    )

    p, n = size(X)

    @inbounds for person in 1:n, snp in 1:p
        if ismissing(X[snp, person])
            #find where snp is located in phase
            hap1_position = searchsortedlast(phase[person].strand1.start, snp)
            hap2_position = searchsortedlast(phase[person].strand2.start, snp)

            #find the correct haplotypes 
            hap1 = phase[person].strand1.haplotypelabel[hap1_position]
            hap2 = phase[person].strand2.haplotypelabel[hap2_position]

            # imputation step 
            X[snp, person] = H[snp, hap1] + H[snp, hap2]
        end
    end

    return nothing
end
