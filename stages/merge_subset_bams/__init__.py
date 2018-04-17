"""
Merges a set of BAMs that are split by locus and barcode into a per-barcode BAM

"""
import martian
import itertools
import os
import subprocess
import tenkit.bam as tk_bam

__MRO__ = """
stage MERGE_BAMS(
    in bam[] subset_bams,
    in string[] subset_bam_keys,
    out bam[] merged_bams,
    src py "merge_bams",
) split using(
    in bam node_bams,
    out bam merged_bam,
)
"""

def split(args):
    chunks = []
    for node_id, bams in itertools.groupby(zip(args.subset_bam_keys, args.subset_bams), key=lambda x: x[0]):
        bams = zip(*list(bams))[1]
        chunks.append({'node_bams': bams, 'node_id': node_id, '__mem_gb': 12*6, '__threads': 12})
    return {'chunks': chunks}


def main(args, outs):
    args.coerce_strings()
    tmp_bam = martian.make_path(str(args.node_id) + '.unsorted.bam')
    tk_bam.concatenate(tmp_bam, args.node_bams)
    outs.merged_bam = martian.make_path('{}.bam'.format(args.node_id))
    subprocess.check_call(['sambamba', 'sort', '-t', str(args.__threads), '-o', outs.merged_bam, tmp_bam])
    os.remove(tmp_bam)


def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    outs.merged_bams = [chunk.merged_bam for chunk in chunk_outs]