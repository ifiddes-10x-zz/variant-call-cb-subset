"""
This stage is a wrapper for GATK
"""
import martian
import os
import re
import subprocess

__MRO__ = '''
stage CALL_SNPINDELS_SUBSET(
    in string targets_file,
    in bam[] subset_bams,
    in string[] barcode_subsets,
    in string[] cell_ids,
    in string[] node_ids,
    in string reference_path,
    in string gatk_path,
    out vcf variants,
    out tsv barcode_map,
    src py "stages/call_snpindels_subset",
) split using(
    in bam subset_bam,
    out vcf subset_variants,
)
'''


def split(args):
    """
    Perform split
    """
    return {'chunks': [{'subset_bam':bam, 'node_id': n, '__mem_gb': 6} for bam, n in zip(args.subset_bams,
                                                                                         args.node_ids)]}


def main(args, outs):
    """
    Wrapper for GATK
    """
    tmp = martian.make_path('{}.vcf'.format(args.node_id))  # node ID in original tree
    ref = os.path.join(args.reference_path, 'fasta', 'genome.fa')
    cmd = ['java', '-jar', args.gatk_path, 'HaplotypeCaller',
           '-R', ref, '-L', args.targets_file,
           '-I', args.subset_bam, '-O', tmp,
           '--native-pair-hmm-threads', '1']
    subprocess.check_call(cmd)
    subprocess.check_call(['bgzip', tmp])
    subprocess.check_call(['tabix', '-p', 'vcf', tmp + '.gz'])
    outs.subset_variants = tmp + '.gz'


def join(args, outs, chunk_defs, chunk_outs):
    """
    Merges the variants produced. For this vcf-merge will be used. Also reports a barcode map.
    """
    outs.coerce_strings()
    vcfs = [x.subset_variants for x in chunk_outs]
    cmd = ['vcf-merge'] + vcfs
    #with open(outs.variants, 'w') as outf:
    tmp = martian.make_path('tmp.vcf')
    with open(tmp, 'w') as outf:
        subprocess.check_call(cmd, stdout=outf)
    # fix the header
    r = re.compile('chnk[0-9A-Za-z-]+\/files\/([0-9]+)_[0-9]+')
    with open(outs.variants, 'w') as outf:
        for l in open(tmp):
            outf.write(r.sub('\\1', l))
    with open(outs.barcode_map, 'w') as outf:
        for node_id, cell_ids, bcodes in sorted(zip(args.node_ids, args.cell_ids, args.barcode_subsets),
                                                key=lambda x: int(x[0])):
            outf.write('\t'.join(map(str, [node_id, cell_ids, bcodes])) + '\n')