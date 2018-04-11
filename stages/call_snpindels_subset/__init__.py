"""
This stage is a wrapper for GATK
"""
import martian
import os
import subprocess

__MRO__ = '''
stage CALL_SNPINDELS_SUBSET(
    in string targets_file,
    in bam[] subset_bams,
    in string[] barcode_subsets,
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
    return {'chunks': [{'subset_bam':bam, 'index': i, '__mem_gb': 6} for i, bam in enumerate(args.subset_bams)]}


def main(args, outs):
    """
    Wrapper for GATK
    """
    tmp = martian.make_path('{}.vcf'.format(args.index))  # use index to map name -> barcode set
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
    with open(outs.variants, 'w') as outf:
        subprocess.check_call(cmd, stdout=outf)
    with open(outs.barcode_map, 'w') as outf:
        for i, bcodes in enumerate(args.barcode_subsets):
            outf.write('\t'.join(map(str, [i, ','.join(bcodes)])) + '\n')