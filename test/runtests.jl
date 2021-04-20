using MendelImpute
using Test
using Random
using JLSO
using VCFTools
using CSV
using Statistics
using BenchmarkTools
using DataFrames
using VariantCallFormat

# change directory to where data is located
cd(normpath(MendelImpute.datadir()))

@testset "search breakpoints" begin
    X1 = Vector{Union{Float64, Missing}}([0,0,0,1,1,1])
    h1 = zeros(6)
    h2 = zeros(6)
    h3 = ones(6)
    bkpt_optim, err_optim = MendelImpute.search_breakpoint(X1, h1, h2, h3)
    @test bkpt_optim == 3
    @test err_optim == 0

    X2 = Vector{Union{Float64, Missing}}([0,0,1,1,1,1])
    bkpt_optim, err_optim = MendelImpute.search_breakpoint(X2, h1, h2, h3)
    @test bkpt_optim == 2
    @test err_optim == 0

    X3 = Vector{Union{Float64, Missing}}([missing,0,1,1,1,1])
    bkpt_optim, err_optim = MendelImpute.search_breakpoint(X3, h1, h2, h3)
    @test bkpt_optim == 2
    @test err_optim == 0

    X4 = Vector{Union{Float64, Missing}}([0,0,0,missing,1,2])
    bkpt_optim, err_optim = MendelImpute.search_breakpoint(X4, h1, h2, h3)
    @test bkpt_optim == 3
    @test err_optim == 1
end

@testset "intersect code" begin
    x = [1, 3, 4, 5, 7, 9]
    y = [2, 3, 5, 6]
    @test MendelImpute.intersect_size_sorted(x, y) == 2
    MendelImpute.intersect_lange!(x, y)
    @test all(x .== [3, 5])
    @test all(y .== [2, 3, 5, 6])

    x = [1, 3, 4, 7]
    y = [2, 3, 5, 6, 7, 10]
    @test MendelImpute.intersect_size_sorted(x, y) == 2
    MendelImpute.intersect_lange!(x, y)
    @test all(x .== [3, 7])
    @test all(y .== [2, 3, 5, 6, 7, 10])

    # allow repeats, although we don't have any in MendelImpute
    x = [3, 4, 7, 7, 7, 10] 
    y = [2, 3, 5, 7, 7, 10]
    @test MendelImpute.intersect_size_sorted(x, y) == 4
    MendelImpute.intersect_lange!(x, y)
    @test all(x .== [3, 7, 7, 10])
    @test all(y .== [2, 3, 5, 7, 7, 10])
end

@testset "dynamic programming" begin
    T = Tuple{Int32, Int32}
    windows = 4
    haplotype_set = [T[] for i in 1:windows]
    push!(haplotype_set[1], (1, 4))
    push!(haplotype_set[1], (1, 5))
    push!(haplotype_set[1], (1, 6))
    push!(haplotype_set[1], (2, 4))
    push!(haplotype_set[1], (2, 5))
    push!(haplotype_set[1], (2, 6))
    push!(haplotype_set[1], (3, 4))
    push!(haplotype_set[1], (3, 5))
    push!(haplotype_set[1], (3, 6))

    push!(haplotype_set[2], (1, 5))
    push!(haplotype_set[2], (1, 7))
    push!(haplotype_set[2], (1, 8))
    push!(haplotype_set[2], (2, 5))
    push!(haplotype_set[2], (2, 7))
    push!(haplotype_set[2], (2, 8))
    push!(haplotype_set[2], (6, 5))
    push!(haplotype_set[2], (6, 7))
    push!(haplotype_set[2], (6, 8))

    push!(haplotype_set[3], (1, 1))
    push!(haplotype_set[3], (1, 5))
    push!(haplotype_set[3], (3, 1))
    push!(haplotype_set[3], (3, 5))

    push!(haplotype_set[4], (4, 4))
    push!(haplotype_set[4], (4, 6))
    push!(haplotype_set[4], (4, 8))
    push!(haplotype_set[4], (5, 4))
    push!(haplotype_set[4], (5, 6))
    push!(haplotype_set[4], (5, 8))
    push!(haplotype_set[4], (8, 4))
    push!(haplotype_set[4], (8, 6))
    push!(haplotype_set[4], (8, 8))

    sol_path, next_pair, subtree_err, best_err = 
        MendelImpute.connect_happairs!(haplotype_set)

    @test best_err == 1.0
    @test sol_path == [(1, 5), (1, 5), (1, 5), (5, 4)]
    @test all(subtree_err[1] .== [2.0; 1.0; 2.0; 3.0; 2.0; 3.0; 5.0; 2.0; 3.0])
    @test all(subtree_err[2] .== [1.0; 2.0; 2.0; 2.0; 5.0; 5.0; 2.0; 5.0; 5.0])
    @test all(subtree_err[3] .== [4.0; 1.0; 4.0; 1.0])
    @test all(subtree_err[end] .== 0.0)
    @test all(next_pair[1] .== [1; 1; 1; 4; 1; 4; 1; 1; 7])
    @test all(next_pair[2] .== 2)
    @test all(next_pair[3] .== [1; 4; 1; 4])
    @test all(next_pair[4] .== 0)
end

@testset "compress" begin
    # perform compression
    reffile = "ref.excludeTarget.vcf.gz"
    tgtfile = "target.typedOnly.masked.vcf.gz"
    outfile = "ref.excludeTarget.jlso"
    @time compress_haplotypes(reffile, tgtfile, outfile)

    # read compressed data
    @time compressed_Hunique = JLSO.load(outfile)[:compressed_Hunique]

    # tests on summary statistics
    @test MendelImpute.nhaplotypes(compressed_Hunique) == 4800
    @test MendelImpute.nwindows(compressed_Hunique) == 16
    @test MendelImpute.max_dim(compressed_Hunique) == (313, 767)
    @test MendelImpute.get_window_widths(compressed_Hunique) == [313, 312, 
        313, 312, 313, 312, 313, 312, 313, 312, 313, 312, 313, 312, 313, 312]
    @test MendelImpute.count_haplotypes_per_window(compressed_Hunique) == [736, 
        752, 733, 742, 765, 650, 666, 745, 679, 695, 682, 767, 720, 665, 711, 
        623]
    @test MendelImpute.avg_haplotypes_per_window(compressed_Hunique) == 708.1875

    # tests on conversion
    @test MendelImpute.unique_idx_to_complete_idx(100, 10, compressed_Hunique) == 127
    @test MendelImpute.unique_all_idx_to_complete_idx(100, 10, compressed_Hunique) == 112
    @test MendelImpute.unique_all_idx_to_unique_typed_idx(100, 10, compressed_Hunique) == 93
    @test MendelImpute.complete_idx_to_unique_all_idx(100, 10, compressed_Hunique) == 91
    @test MendelImpute.complete_idx_to_unique_typed_idx(100, 10, compressed_Hunique) == 85
end

@testset "impute & phasing" begin
    Random.seed!(2020)
    tgtfile = "target.typedOnly.masked.vcf.gz"
    reffile = "ref.excludeTarget.jlso"
    outfile = "imputed.vcf.gz"
    @time phase(tgtfile, reffile, outfile)

    # check error rate
    Ximputed = convert_gt(Float64, "imputed.vcf.gz")
    Xtrue = convert_gt(Float64, "target.full.vcf.gz")
    m, n = size(Xtrue)
    error_rate = sum(Xtrue .!= Ximputed) / m / n
    @test (m, n) == (100, 36063)
    @test error_rate ≈ 0.0007129190583146162

    # check per-sample error
    quality = CSV.read("imputed.sample.error", DataFrame)
    @test size(quality, 1) == 100
    @test count(isone, quality[!, :error]) == 65
    @test mean(quality[!, :error]) ≈ 0.9996911073712372
    @test std(quality[!, :error]) ≈ 0.0006256546642253605

    # check per-SNP error
    reader = VCF.Reader(openvcf(outfile, "r"))
    snpscores = Vector{Float64}(undef, nrecords(outfile))
    for (i, record) in enumerate(reader)
        snpscores[i] = parse(Float64, VCF.info(record)[1].second)
    end
    close(reader)
    @test length(snpscores) == 36063
    @test count(isone, snpscores) == 29934
    @test mean(snpscores) ≈ 0.9997004408951002
    @test std(snpscores) ≈ 0.0007984939339310984
end

@testset "PLINK and dosage inputs" begin
    # compute .vcf.gz (hard genotype) error rate
    tgtfile = "target.typedOnly.masked.vcf.gz"
    reffile = "ref.excludeTarget.jlso"
    outfile = "imputed.vcf.gz"
    @time phase(tgtfile, reffile, outfile);
    Ximputed = convert_gt(Float64, "imputed.vcf.gz")
    Xtrue = convert_gt(Float64, "target.full.vcf.gz")
    m, n = size(Xtrue)
    vcf_error = sum(Xtrue .!= Ximputed) / m / n

    # PLINK error rate
    tgtfile = "target.typedOnly.masked" # do no include bim/bed/fam in filename
    reffile = "ref.excludeTarget.jlso"
    outfile = "imputed.vcf.gz"
    @time phase(tgtfile, reffile, outfile)
    Ximputed = convert_gt(Float64, "imputed.vcf.gz")
    plink_error = sum(Xtrue .!= Ximputed) / m / n
    @test plink_error ≈ vcf_error

    # dosage error rate
    tgtfile = "target.typedOnly.dosages.masked.vcf.gz"
    reffile = "ref.excludeTarget.jlso"
    outfile = "imputed.vcf.gz"
    @time phase(tgtfile, reffile, outfile, dosage=true)
    Ximputed = convert_gt(Float64, "imputed.vcf.gz")
    dosage_error = sum(Xtrue .!= Ximputed) / m / n
    @test round(dosage_error, digits = 4) == round(vcf_error, digits=4)
end

@testset "decompress ultra-compressed files" begin
    # saving
    tgtfile = "target.typedOnly.masked.vcf.gz"
    reffile = "ref.excludeTarget.jlso"
    outfile = "imputed.jlso" # outfile ends in .jlso
    @time phase(tgtfile, reffile, outfile);

    # loading requires original reference panel in vcf format
    tgtfile = "imputed.jlso"
    reffile = "ref.excludeTarget.vcf.gz"
    X1, X2, phaseinfo, sampleID, H = convert_compressed(Float64, 
        tgtfile, reffile);

    # test error rate
    Xtrue = convert_gt(Float64, "target.full.vcf.gz")
    Ximputed = (X1 + X2)' # transpose X1 and X2
    m, n = size(Ximputed)
    error_rate = sum(Xtrue .!= Ximputed) / m / n
    @test error_rate ≈ 0.0007129190583146162
end

@testset "0 allocations in phasing" begin
    reffile = "ref.excludeTarget.jlso"
    compressed_Hunique = JLSO.load(reffile)[:compressed_Hunique]

    # first person's optimal haplotype in each window (complete index)
    happair1_original = Int32[56, 624, 22, 20, 54, 221, 515, 172, 989, 19, 
        171, 171, 570, 37, 570, 106]
    happair2_original = Int32[624, 989, 624, 624, 624, 515, 989, 989, 4115,
        171, 989, 202, 1681, 1681, 4115, 570]

    # preallocate vectors
    survivors1=Int32[]
    survivors2=Int32[]
    sizehint!(survivors1, 60000)
    sizehint!(survivors2, 60000)
    happair1 = copy(happair1_original)
    happair2 = copy(happair2_original)
    MendelImpute.phase_sample!(happair1, happair2,
        compressed_Hunique, survivors1, survivors2)

    # zero allocation and memory
    b = @benchmark MendelImpute.phase_sample!(happair1, happair2, 
        $compressed_Hunique, $survivors1, $survivors2) setup=(happair1=
        copy($happair1_original);happair2=copy($happair2_original))

    @test b.allocs == 0
    @test b.memory == 0
end

@testset "0 allocation in haplopair" begin
    tgtfile = "target.typedOnly.masked.vcf.gz"
    reffile = "ref.excludeTarget.jlso"
    compressed_Hunique = JLSO.load(reffile)[:compressed_Hunique]

    # a bunch of setup code
    X, X_sampleID, X_chr, X_pos, X_ids, X_ref, X_alt = VCFTools.convert_gt(
        UInt8, tgtfile, trans=true, save_snp_info=true, msg = 
        "Importing genotype file...");
    people = size(X, 2)
    tgt_snps = size(X, 1)
    ref_snps = length(compressed_Hunique.pos)
    threads = Threads.nthreads()
    max_width, max_d = MendelImpute.max_dim(compressed_Hunique)
    inv_sqrt_allele_var = nothing
    overlap = compressed_Hunique.overlap
    windows = MendelImpute.nwindows(compressed_Hunique)
    num_unique_haps = round(Int, 
        MendelImpute.avg_haplotypes_per_window(compressed_Hunique))

    # preallocated arrays in phase
    ph = [MendelImpute.HaplotypeMosaicPair(ref_snps) for i in 1:people]
    haplotype1 = [zeros(Int32, windows) for i in 1:people]
    haplotype2 = [zeros(Int32, windows) for i in 1:people]
    
    # working arrys in happair
    happair1 = [ones(Int32, people)               for _ in 1:threads]
    happair2 = [ones(Int32, people)               for _ in 1:threads]
    hapscore = [zeros(Float32, people)            for _ in 1:threads]
    Xwork    = [zeros(Float32, max_width, people) for _ in 1:threads]
    Hwork    = [zeros(Float32, max_width, max_d)  for _ in 1:threads]

    # window 1
    w = 1
    Hw_aligned = compressed_Hunique.CW_typed[w].uniqueH
    Xrange = MendelImpute.extend_to_overlap_range(compressed_Hunique, w, overlap)
    Xw_aligned = view(X, Xrange, :)
    d = size(Hw_aligned, 2)
    id = Threads.threadid()

    # global search must allocate M, Hwork, N for each window
    # but there's no other allocation
    M = zeros(Float32, size(Hw_aligned, 2), size(Hw_aligned, 2))
    N = zeros(Float32, size(Xw_aligned, 2), size(Hw_aligned, 2))
    Hwork = convert(Matrix{Float32}, Hw_aligned)
    b = @benchmark MendelImpute.haplopair!($Xw_aligned, $Hw_aligned, 
        happair1=$(happair1[id]), happair2=$(happair2[id]), 
        hapscore=$(hapscore[id]), Xwork=$(Xwork[id]), M=$M, N=$N, Hwork=$Hwork)

    @test b.allocs == 0
    @test b.memory == 0
end
