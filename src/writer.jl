function Base.write(
	writer::VCF.Writer,
	phase::Vector{HaplotypeMosaicPair},
    H::AbstractMatrix,
	)
	n = 0
	samples = length(phase)
	records = phase[1].strand1.length
	for i in 1:records
		for j in 1:samples
            #find where snp is located in phase
            # println("i = $i, j = $j")
            hap1_position = searchsortedlast(phase[j].strand1.start, i)
            hap2_position = searchsortedlast(phase[j].strand2.start, i)

            #find the correct haplotypes 
            hap1 = phase[j].strand1.haplotypelabel[hap1_position]
            hap2 = phase[j].strand2.haplotypelabel[hap2_position]

            # actual allele
			a1 = H[j, hap1]
			a2 = H[j, hap2]

			# write to io
			n += write(writer.stream, a1, '|', a2, '\t')
		end
		n += write(writer.stream, '\n')
	end
	return n
end
