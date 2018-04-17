"""
This stage is a wrapper for GATK
"""
import martian
import os
import collections
import subprocess
import tenkit.tabix as tk_tabix
from vclib.locus import build_loci
__MRO__ = '''
stage CALL_SNPINDELS_SUBSET(
    in string targets_file,
    in bam[] merged_bams,
    in string[] barcode_subsets,
    in string[] node_ids,
    in string reference_path,
    in string gatk_path,
    out vcf variants,
    out tsv barcode_map,
    src py "stages/call_snpindels_subset",
) split using(
    in bam subset_bam,
    in string locus,
    out vcf subset_variants,
)
'''


def split(args):
    """
    Perform split. If targets_file is set, then no locus split will be applied.
    """
    if args.targets_file is not None:
        return {'chunks': [{'subset_bam': bam, 'node_id': n, '__mem_gb': 6} for bam, n in zip(args.merged_bams,
                                                                                             args.node_ids)]}
    else:
        loci = build_loci(args.reference_path)
        chunks = []
        for bam, node_id in zip(args.merged_bams, args.node_ids):
            for locus in loci:
                chunks.append({'locus': locus, 'subset_bam': bam, 'node_id': node_id,
                               '__mem_gb': 6 * 6, '__threads': 6})
        return {'chunks': chunks, 'join': {'__mem_gb': 16}}


def main(args, outs):
    """
    Wrapper for GATK
    """
    tmp = martian.make_path('tmp.vcf')
    ref = os.path.join(args.reference_path, 'fasta', 'genome.fa')
    cmd = ['java', '-jar', args.gatk_path, 'HaplotypeCaller',
           '-R', ref,
           '-I', args.subset_bam, '-O', tmp,
           '--native-pair-hmm-threads', str(args.__threads)]
    if args.locus is None:
        assert args.targets_file is not None
        cmd.extend(['-L', args.targets_file])
    else:
        bed = martian.make_path('region.bed')
        with open(bed, 'w') as outf:
            outf.write(args.locus)
        cmd.extend(['-L', bed])
    subprocess.check_call(cmd)

    # fix the name
    tmp2 = martian.make_path('{}.vcf'.format(args.node_id))  # node ID in original tree
    with open(tmp2, 'w') as outf:
        for l in open(tmp):
            if l.startswith('#CHROM'):
                l = l.split()
                l[-1] = args.node_id
                l = '\t'.join(l) + '\n'
            outf.write(l)

    subprocess.check_call(['bgzip', tmp2])
    subprocess.check_call(['tabix', '-p', 'vcf', tmp2 + '.gz'])
    outs.subset_variants = tmp2 + '.gz'


def join(args, outs, chunk_defs, chunk_outs):
    """
    Merges the variants produced. For this vcf-merge will be used. Also reports a barcode map.
    If no targets file was produced, first concatenates VCFs across the loci.
    """
    outs.coerce_strings()

    if args.targets_file is None:
        # concatenate VCFs
        vcfs_by_node = collections.defaultdict(list)
        for o, d in zip(chunk_outs, chunk_defs):
            vcfs_by_node[d.node_id].append(o.subset_variants)

        cat_vcfs = {}
        for node_id, vcfs in vcfs_by_node.iteritems():
            cmd = ['vcf-concat'] + vcfs
            tmp_vcf = martian.make_path('{}.vcf'.format(node_id))
            sorted_vcf = martian.make_path('{}.sorted.vcf'.format(node_id))
            with open(tmp_vcf, 'w') as outf:
                subprocess.check_call(cmd, stdout=outf)
            tk_tabix.sort_vcf(tmp_vcf, sorted_vcf)
            tk_tabix.index_vcf(sorted_vcf)
            assert os.path.exists(sorted_vcf + '.gz')
            cat_vcfs[node_id] = sorted_vcf + '.gz'

        vcfs = [cat_vcfs[x] for x in args.node_ids]
    else:
        vcfs = [x.subset_variants for x in chunk_outs]

    if len(vcfs) > 1:
        cmd = ['vcf-merge'] + vcfs
        with open(outs.variants, 'w') as outf:
            subprocess.check_call(cmd, stdout=outf)
    else:
        os.rename(vcfs[0], outs.variants)

    with open(outs.barcode_map, 'w') as outf:
        for node_id, bcodes in sorted(zip(args.node_ids, args.barcode_subsets), key=lambda x: int(x[0])):
            outf.write('\t'.join(map(str, [node_id, bcodes])) + '\n')
