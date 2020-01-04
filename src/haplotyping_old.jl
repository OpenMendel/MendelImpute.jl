"""
    phase(X, H, width=400, verbose=true)

Phasing (haplotying) of genotype matrix `X` from a pool of haplotypes `H`
by sliding windows. This is Hua's code. 

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
