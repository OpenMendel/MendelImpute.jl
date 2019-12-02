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
    H::AbstractMatrix{T};
    width::Int    = 700,
    verbose::Bool = true
    ) where T <: Real

    # declare some constants
    snps, people = size(X)
    haplotypes = size(H, 2)
    windows = floor(Int, snps / width)

    # get redundant haplotype sets. 
    hapset = compute_optimal_halotype_set(X, H, width=width, verbose=verbose)

    # allocate working arrays
    phase = [HaplotypeMosaicPair(snps) for i in 1:people]
    haplo_chain = ([copy(hapset[i].strand1[1]) for i in 1:people], [copy(hapset[1].strand2[1]) for i in 1:people])
    chain_next  = (BitVector(undef, haplotypes), BitVector(undef, haplotypes))
    window_span = (ones(Int, people), ones(Int, people))

    # TODO: parallel computing
    # begin intersecting haplotypes window by window 
    @inbounds for i in 1:people, w in 2:windows

        # decide whether to cross over based on the larger intersection
        chain_next[1] .= haplo_chain[1][i] .& hapset[i].strand1[w] # not crossing over
        chain_next[2] .= haplo_chain[1][i] .& hapset[i].strand2[w] # crossing over
        if sum(chain_next[1]) < sum(chain_next[2])
            hapset[i].strand1[w], hapset[i].strand2[w] = hapset[i].strand2[w], hapset[i].strand1[w]
        end        

        # strand 1 
        chain_next[1] .= haplo_chain[1][i] .& hapset[i].strand1[w]
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

        # strand 2
        chain_next[2] .= haplo_chain[2][i] .& hapset[i].strand2[w]
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

    # handle last few windows separately, since intersection may not become empty
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
        push!(phase[i].strand1.start, 1)
        push!(phase[i].strand1.haplotypelabel, hap1)
        push!(phase[i].strand2.start, 1)
        push!(phase[i].strand2.haplotypelabel, hap2)
    end

    # find optimal break points and record info to phase. 
    # TODO: handle last window separately since view() on X or H is not complete
    strand1_intersect = chain_next[1]
    strand2_intersect = chain_next[2]
    for i in 1:people, w in 2:windows
        
        strand1_intersect .= hapset[i].strand1[w - 1] .& hapset[i].strand1[w]
        if sum(strand1_intersect) == 0
            # search breakpoints
            Xi = view(X, ((w - 2) * width + 1):(w * width), i)
            Hi = view(H, ((w - 2) * width + 1):(w * width), :)
            s2 = findfirst(hapset[i].strand2[w]) :: Int64
            s1_cur  = findfirst(hapset[i].strand1[w - 1]) :: Int64
            s1_next = findfirst(hapset[i].strand1[w]) :: Int64
            bkpt, err_optim = search_breakpoint(Xi, Hi, s2, (s1_cur, s1_next))

            # record info into phase
            push!(phase[i].strand1.start, (w - 2) * width + 1 + bkpt)
            push!(phase[i].strand1.haplotypelabel, s1_next)
        end

        strand2_intersect .= hapset[i].strand2[w - 1] .& hapset[i].strand2[w]
        if sum(strand2_intersect) == 0
            # search breakpoints
            Xi = view(X, ((w - 2) * width + 1):(w * width), i)
            Hi = view(H, ((w - 2) * width + 1):(w * width), :)
            s1 = findfirst(hapset[i].strand1[w]) :: Int64
            s2_cur  = findfirst(hapset[i].strand2[w - 1]) :: Int64
            s2_next = findfirst(hapset[i].strand2[w]) :: Int64
            bkpt, err_optim = search_breakpoint(Xi, Hi, s1, (s2_cur, s2_next))

            # record info into phase
            push!(phase[i].strand2.start, (w - 2) * width + 1 + bkpt)
            push!(phase[i].strand2.haplotypelabel, s2_next)
        end
    end

    return hapset, phase 
end

"""
    impute!(X, H, phase)

Imputes `X` completely using segments of haplotypes `H` where segments are stored in `phase`. 
Non-missing entries in `X` can be different after imputation. 
"""
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

"""
    impute2!(X, H, phase)

Imputes missing entries of `X` using corresponding haplotypes `H` via `phase` information. 
Non-missing entries in `X` will not change. 
"""
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
