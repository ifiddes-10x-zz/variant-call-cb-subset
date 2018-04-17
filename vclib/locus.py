import longranger.cnv.contig_manager as contig_manager


def build_loci(ref_path, chunk_size=5 * 10 ** 7):
    """
    Given a ref path, builds a set of loci
    :param ref_path: Path to cellranger reference
    :param chunk_size: Number of bases to pack into one bin
    :return:
    """
    ref = contig_manager.contig_manager(ref_path)
    # greedy implementation of the bin packing problem
    loci = []
    this_locus = []
    bin_size = 0
    for chrom, chrom_length in ref.contig_lengths.iteritems():
        region_start = 0
        while region_start < chrom_length:
            start = region_start
            end = min(region_start + chunk_size, chrom_length)
            this_locus.append([chrom, start, end])
            bin_size += end - start
            if bin_size >= chunk_size:
                loci.append('\n'.join(['\t'.join(map(str, x)) for x in this_locus]) + '\n')
                this_locus = []
                bin_size = 0
            region_start = end
    return loci
