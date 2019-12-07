"""
    write(writer, phase, H)

Writes phased data to file.

TODO: Make this efficient (current takes 2x time than phasing and 100x memory)
TODO: Check if `a1` and `a2` needs to change depending on REF/ALT. 
"""
function Base.write(
    writer::VCF.Writer,
    phase::Vector{HaplotypeMosaicPair},
    H::AbstractMatrix,
    )
    samples = length(phase)
    records = phase[1].strand1.length
    for i in 1:records
        for j in 1:samples
            #find where snp is located in phase
            hap1_position = searchsortedlast(phase[j].strand1.start, i)
            hap2_position = searchsortedlast(phase[j].strand2.start, i)

            #find the correct haplotypes 
            hap1 = phase[j].strand1.haplotypelabel[hap1_position]
            hap2 = phase[j].strand2.haplotypelabel[hap2_position]

            # actual allele
            a1 = convert(Int, H[i, hap1])
            a2 = convert(Int, H[i, hap2])

            # write to io
            print(writer.stream, a1, '|', a2, '\t')
        end
        write(writer.stream, '\n')
    end
    flush(writer)
end
