"""
    make_refvcf_file(H, filename, phased)

Creates a .vcf file with `filename` based on reference panels `H`. Consecutive columns are treated as genotypes.
REF/ALT alleles are always A/C.

# Inputs:
+ `H`: BitMatrix of haplotypes. Each column is a haplotype. e.g. Columns 1 and 2 form the genotype for sample 1. 
+ `filename`: A string for the resulting .vcf file. 
+ `phased`: True uses '|' as separator. False uses '/' as separator. 
"""
function make_refvcf_file(
    H::AbstractMatrix;
    vcffilename::AbstractString = "simulated_ref.vcf", 
    phased::Bool = true,
    marker_chrom::Vector{String} = ["1" for i in 1:size(H, 1)],
    marker_pos::Vector{Int} = collect(1:size(H, 1)),
    marker_ID::Vector{String} = ["tgt_snp_$i" for i in 1:size(H, 1)],
    marker_REF::Vector{String} = ["A" for i in 1:size(H, 1)],
    marker_ALT::Vector{String} = ["C" for i in 1:size(H, 1)]
    )

    p, d = size(H)
    separator = (phased ? '|' : '/')
    iseven(d) || error("make_vcf_file: number of haplotypes must be even but was $d")

    # first write minimal meta information
    io = openvcf(vcffilename, "w")
    write(io, "##fileformat=VCFv4.3\n")
    write(io, "##source=MendelImpute\n")
    write(io, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

    # header line
    write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    for i in 1:Int(d / 2)
        write(io, "\tref$i")
    end
    write(io, "\n")

    # write phase info
    for snp in 1:p
        write(io, "1\t$snp\tref_snp_$snp\tA\tC\t.\tPASS\t.\tGT")
        for i in 1:2:d
            write(io, string("\t", Int(H[snp, i]), separator, Int(H[snp, i + 1])))
        end
        write(io, "\n")
    end
    close(io)
end

"""
    make_tgtvcf_file(X, filename)

Creates a .vcf file given a genotype matrix `X`. Missing entries in `X` will become `./.`. 
Otherwise, X[i, j] should be either 0, 1, or 2. 0 = '0/0', 1 = '1/0', and 2 = '1/1'. REF/ALT 
alleles are always A/C.

# Inputs:
+ `X`: Matrix of 0, 1, or 2. Each column is a person's genotype. 
+ `filename`: A string for the resulting .vcf file. 
"""
function make_tgtvcf_file(
    X::AbstractMatrix;
    vcffilename::AbstractString = "simulated_tgt.vcf", 
    phased::Bool = false,
    marker_chrom::Vector{String} = ["1" for i in 1:size(X, 1)],
    marker_pos::Vector{Int} = collect(1:size(X, 1)),
    marker_ID::Vector{String} = ["tgt_snp_$i" for i in 1:size(X, 1)],
    marker_REF::Vector{String} = ["A" for i in 1:size(X, 1)],
    marker_ALT::Vector{String} = ["C" for i in 1:size(X, 1)]
    )

    p, d = size(X)
    lc = length(marker_chrom)
    lp = length(marker_pos)
    li = length(marker_ID)
    lr = length(marker_REF)
    la = length(marker_ALT)
    p == lc == lp == li == lr == la || error("There are $p markers in X but CHROM/POS/ID/REF/ALT vectors are of length $lc, $lp, $li, $lr, $la")
    separator = (phased ? '|' : '/')

    # first write minimal meta information
    io = openvcf(vcffilename, "w")    
    write(io, "##fileformat=VCFv4.3\n")
    write(io, "##source=MendelImpute\n")
    write(io, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

    # header line
    write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    for i in 1:d
        write(io, "\ttarget$i")
    end
    write(io, "\n")

    # write phase info
    for snp in 1:p
        write(io, marker_chrom[snp], '\t', string(marker_pos[snp]), '\t', marker_ID[snp], '\t', 
            marker_REF[snp], '\t', marker_ALT[snp], "\t.\tPASS\t.\tGT")
        @inbounds for i in 1:d
            if ismissing(X[snp, i])
                genotype = "./."
            elseif X[snp, i] == 2
                genotype = "1/1"
            elseif X[snp, i] == 1
                genotype = "1/0"
            elseif X[snp, i] == 0
                genotype = "0/0"
            else
                error("genotypes can only be 0, 1, 2, or missing, but was $(X[snp, i])")
            end
            write(io, '\t', genotype)
        end
        write(io, '\n')
    end
    close(io)
end

"""
    simulate_markov_haplotypes(p, d, prob)

Simulates a haplotype matrix as a markov chain. The `i`th allele (0 or 1) will transition 
to the opposite allele (0 or 1) with probability `prob` at the `i + 1`th allele. 

# Inputs
- `p`: Length of each haplotype.
- `d`: Total number of haplotypes. 
- `prob`: transition probability.

# Output
- `H`: `p x d` haplotype matrix. Each column is a haplotype
"""
function simulate_markov_haplotypes(
    p::Int64, 
    d::Int64;
    prob::AbstractFloat = 0.25,
    vcffilename::String = ""
    )
    @assert 0 < prob < 1 "transition probably `prob` should be between 0 and 1, got $prob"

    H = falses(p, d)
    @inbounds for j in 1:d
        H[1, j] = rand(Bool)
        for i in 2:p
            H[i, j] = (rand() < prob ? !H[i - 1, j] : H[i - 1, j])
        end
    end

    # create VCF file on disk if vcffilename is provided
    if vcffilename != ""
        make_refvcf_file(H, vcffilename=vcffilename)
    end

    return H
end

"""
    simulate_uniform_haplotypes(p, d, prob)

Simulates a haplotype matrix `H` where `H[i, j] = 1` with probability `prob`. 

# Inputs
- `p`: Length of each haplotype.
- `d`: Total number of haplotypes. 
- `prob`: probability that an entry in H is 1.

# Output
- `H`: `p x d` haplotype matrix. Each column is a haplotype
"""
function simulate_uniform_haplotypes(
    p::Int64, 
    d::Int64;
    prob = 0.25,
    )
    @assert 0 < prob < 1 "prob should be between 0 and 1, got $prob"

    H = falses(p, d)
    @inbounds for j in 1:d, i in 1:p
        if rand() < prob
            H[i, j] = true
        end
    end
    return H
end

"""
    simulate_genotypes(H; block_length)

Simulates a genotype matrix `X` from a pool of haplotypes `H`. Each person's
genotype are divided into contiguous segments, and 2 haplotypes are randomly
chosen from a pool of haplotypes `H` to form the genotype in that segment. 

# Arguments:
- `H`: `p x d` haplotype matrix. Each column is a haplotype. 
- `people`: number of samples
- `T`: Type of output matrix. 
- `min_cross_over`: Minimum number of breakpoints for each person's genotype. 
- `max_cross_over`: Maximum number of breakpoints for each person's genotype. 

# Output:
* `X`: `p x people` genotype matrix. Each column is a person's genotype. 
"""
function simulate_genotypes(
    H::AbstractMatrix,
    people::Int;
    T::Type = Int,
    min_cross_over::Int64=1,
    max_cross_over::Int64=5,
    )
    
    p, d = size(H)
    X = zeros(Union{T, Missing}, p, people)
    min_cross_over <= max_cross_over || error("Please supply min_cross_over and max_cross_over satisfying min_cross_over <= max_cross_over.")

    # loop through each person
    segments = UnitRange{Int64}[]
    sizehint!(segments, max_cross_over)
    for i in 1:people
        # simulate cross overs. no crossovers in first/last windows and 2 crossover cannot occur within 1 window
        cross_overs = rand(min_cross_over:max_cross_over)
        cross_over_location = collect(1:cross_overs)
        while true
            cross_over_location .= sample((width + 1):(p - width - 1), cross_overs, replace=false) 
            sort!(cross_over_location)
            cross_overs == 1 && break
            minimum(diff(cross_over_location)) > width && break
        end
        #create various segments vased on cross over points
        empty!(segments)
        push!(segments, 1:cross_over_location[1])
        for j in 1:(length(cross_over_location) - 1)
            push!(segments, (cross_over_location[j] + 1):cross_over_location[j + 1])
        end
        push!(segments, (cross_over_location[end] + 1):p)
        # fill X with sum of 2 randomly chosen haplotypes in each segment
        for cur_range in segments
            h1, h2 = rand(1:d), rand(1:d)
            X[cur_range, i] .= convert.(T, H[cur_range, h1] .+ H[cur_range, h2])
        end
    end
    return X
end

"""
Simulates genotype from haplotype reference panels. 

# Inputs
+ `H`: Haplotype matrix. Columns are haplotypes
+ `people`: Integer for number of samples you want to simulate. 

Returns `hap_mosaics` and `hap_mosaic_range` where 
+ `hap_mosaics[i]` is a vector haplotypes [(hi, hj), ...] that was used to simulate person `i`'s genotype
+ `hap_mosaic_range` is a vector of ranges (e.g. [1:100, ...]) that records the range of SNPs where (hi, hj) filled.
"""
function simulate_phased_genotypes(
    H::Union{AbstractMatrix, String},
    people::Int;
    T::Type = Int,
    min_cross_over::Int64=1,
    max_cross_over::Int64=5,
    vcffilename::String = "",
    width=400
    )
    if typeof(H) == String
        H = convert_ht(Float32, H)
    end

    p, d = size(H)
    separator = '|'
    iseven(d) || error("simulate_phased_genotypes: number of haplotypes must be even but was $d")

    # create target matrix, where 2 columns form 1 genotype
    X = zeros(Union{T, Missing}, p, 2people)
    hap_mosaics = [Tuple{Int, Int}[] for i in 1:people] 
    hap_mosaic_range = [UnitRange{Int64}[] for i in 1:people]

    segments = UnitRange{Int64}[]
    sizehint!(segments, max_cross_over)
    for i in 1:people
        # simulate cross overs. no crossovers in first/last windows and 2 crossover cannot occur within 1 window
        cross_overs = rand(min_cross_over:max_cross_over)
        cross_over_location = collect(1:cross_overs)
        while true
            cross_over_location .= sample((width + 1):(p - width - 1), cross_overs, replace=false) 
            sort!(cross_over_location)
            cross_overs == 1 && break
            minimum(diff(cross_over_location)) > width && break
        end
        #create various segments vased on cross over points
        empty!(segments)
        push!(segments, 1:cross_over_location[1])
        for j in 1:(length(cross_over_location) - 1)
            push!(segments, (cross_over_location[j] + 1):cross_over_location[j + 1])
        end
        push!(segments, (cross_over_location[end] + 1):p)
        # fill X with 2 randomly chosen haplotypes in each segment
        for cur_range in segments
            h1, h2 = rand(1:d), rand(1:d)
            X[cur_range, 2i - 1] .= convert.(T, H[cur_range, h1])
            X[cur_range, 2i    ] .= convert.(T, H[cur_range, h2])

            # record haplotypes and range
            push!(hap_mosaics[i], (h1, h2))
            push!(hap_mosaic_range[i], cur_range)
        end
    end

    if vcffilename != ""
        make_refvcf_file(X, vcffilename=vcffilename)
    end

    return X, hap_mosaics, hap_mosaic_range
end

function unphase(
    tgtfile::AbstractString;
    outfile::AbstractString = "unphase." * tgtfile
    )

    # create VCF reader and writer
    reader = VCF.Reader(openvcf(tgtfile, "r"))
    writer = VCF.Writer(openvcf(outfile, "w"), header(reader))
    pmeter = Progress(nrecords(tgtfile), 1, "Creating $outfile...")

    # loop over each record (snp)
    for (i, record) in enumerate(reader)
        gtkey = VCF.findgenokey(record, "GT")
        if !isnothing(gtkey) 
            # loop over samples
            for (j, geno) in enumerate(record.genotype)
                # unphased data has separator '/'
                record.data[geno[gtkey][2]] = 0x2f 

                # change heterozygotes to 1/0
                a1 = record.data[geno[gtkey][1]]
                a2 = record.data[geno[gtkey][3]]
                if (a1 == 0x30 && a2 != 0x30) || (a1 != 0x30 && a2 == 0x30)
                    #"0" (REF) => 0x30, "1" (ALT) => 0x31
                    record.data[geno[gtkey][1]] = 0x31
                    record.data[geno[gtkey][3]] = 0x30
                end
            end
        end
        write(writer, record)
        next!(pmeter) #update progress
    end

    # close 
    flush(writer); close(reader); close(writer)
    return nothing
end

