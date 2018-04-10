"""
This stage takes a possorted BAM, a list of barcode subsets produced by PARSE_INPUTS,
"""
import pysam
import tenkit.bam as tk_bam

__MRO__ = '''
stage FILTER_BAM(
    in bam possorted_bam,
    in string[] barcode_subsets,
    out bam[] subset_bams,
    src py "stages/filter_bam_by_cb",
) split using(
    in string barcode_subset,
    out bam subset_bam,
    out bam.bai subset_bam_index,
)
'''


def split(args):
    """
    Loads the barcode subsets produced by PARSE_INPUTS and creates a subset job for each set of barcodes
    """
    # TODO: this will be slow unless I do the below. I don't see an easy way to do this, however.
    # TODO: I would need to produce nested stages somehow.
    #bam_in = tk_bam.create_bam_infile(args.possorted_bam)
    #chunks = tk_bam.chunk_bam_records(bam_in, chunk_bound_key=None, chunk_size_gb=8.0)
    #for c in chunks:
    #    c['__mem_gb'] = 3
    #
    #return {'chunks': chunks, 'join': {'__mem_gb': 6*8, '__threads': 8}}
    return {'chunks': [{'barcode_subset': x} for x in args.barcode_subsets]}


def main(args, outs):
    """
    Given a set of barcodes and a possorted bam, return a new BAM that only contains those barcodes
    """
    useful_bcs = set(args.barcode_subset)

    bam_h = pysam.Samfile(args.possorted_bam)
    outf_h = pysam.Samfile(outs.subset_bam, 'wb', template=bam_h)
    for rec in bam_h:
        try:
            cb = rec.get_tag('CB')
        except KeyError:
            continue
        if cb in useful_bcs:
            outf_h.write(rec)
    outf_h.close()
    tk_bam.index(outs.subset_bam)


def join(args, outs, chunk_defs, chunk_outs):
    """
    Returns the set of BAMs produced
    """
    outs.coerce_strings()
    outs.subset_bams = [x.subset_bam for x in chunk_outs]