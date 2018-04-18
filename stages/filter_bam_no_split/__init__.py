"""
This stage takes a possorted BAM, a list of barcode subsets produced by PARSE_INPUTS and generates BAMs
without splitting across loci. Removes the need for MERGE_SUBSET_BAMS.
"""
import tenkit.bam as tk_bam

__MRO__ = '''
stage FILTER_BAM_NO_SPLIT(
    in bam possorted_bam,
    in string reference_path,
    in string[] barcode_subsets,
    out bam[] subset_bams,
    src py "stages/filter_bam_no_split",
) split using(
    in string cluster_bcs,
    out bam subset_bam,
)
'''


def split(args):
    """
    Loads the barcode subsets produced by PARSE_INPUTS and creates a subset job for each barcode subset
    """
    bc_chunks = []
    for bc_subset in args.barcode_subsets:
            bc_chunks.append({'cluster_bcs': bc_subset, '__mem_gb': 6})
    return {'chunks': bc_chunks}


def main(args, outs):
    """
    Given a set of barcodes and a bam chunk, return a new BAM that only contains those barcodes
    """
    useful_bcs = set(args.cluster_bcs.split(','))

    in_bam = tk_bam.create_bam_infile(args.possorted_bam)
    out_bam, _ = tk_bam.create_bam_outfile(outs.subset_bam, None, None, template=in_bam)
    for rec in in_bam:
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
    outs.subset_bams = [chunk.subset_bam for chunk in chunk_outs]