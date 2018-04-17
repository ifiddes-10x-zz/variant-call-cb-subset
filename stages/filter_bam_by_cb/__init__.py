"""
This stage takes a possorted BAM, a list of barcode subsets produced by PARSE_INPUTS,
"""
import tenkit.bam as tk_bam
from vclib.locus import build_loci

__MRO__ = '''
stage FILTER_BAM(
    in bam possorted_bam,
    in string reference_path,
    in string[] barcode_subsets,
    in string[] node_ids,
    out bam[] subset_bams,
    out string[] subset_bam_keys,
    src py "stages/filter_bam_by_cb",
) split using(
    in string cluster_bcs,
    in string[] locus,
    out bam subset_bam,
)
'''


def split(args):
    """
    Loads the barcode subsets produced by PARSE_INPUTS and creates a subset job for each pair of locus-barcode subset
    """
    # construct loci
    loci = build_loci(args.reference_path)
    # nest loci with clusters
    bc_chunks = []
    for bc_subset, node_id in zip(args.barcode_subsets, args.node_ids):
        for locus in loci:
            bc_chunks.append({'cluster_bcs': bc_subset, 'node_id': node_id, 'locus': locus, '__mem_gb': 6})
    return {'chunks': bc_chunks}


def main(args, outs):
    """
    Given a set of barcodes and a bam chunk, return a new BAM that only contains those barcodes
    """
    useful_bcs = set(args.cluster_bcs.split(','))

    in_bam = tk_bam.create_bam_infile(args.possorted_bam)
    out_bam, _ = tk_bam.create_bam_outfile(outs.subset_bam, None, None, template=in_bam)

    locus = [x.split('\t') for x in args.locus.split('\n') if x]  # strip trailing newline
    for chrom, start, stop in locus:
        for rec in in_bam.fetch(str(chrom), int(start), int(stop)):
            try:
                cb = rec.get_tag('CB')
            except KeyError:
                continue
            if cb in useful_bcs:
                out_bam.write(rec)
    out_bam.close()
    tk_bam.index(outs.subset_bam)


def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()

    outs.subset_bam_keys = []
    outs.subset_bams = []
    for d, o in zip(chunk_defs, chunk_outs):
        outs.subset_bam_keys.append(d.node_id)
        outs.subset_bams.append(o.subset_bam)
