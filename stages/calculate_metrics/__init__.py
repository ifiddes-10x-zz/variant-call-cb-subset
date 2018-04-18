"""
Calculates metrics across the BAM files

metrics:
1. Average depth (meanI)
2. Variance (varI)
3. DPCV (CVI)
4. Breadth of coverage (Number of bases in genome with a read covering)
"""
import collections
import json
import subprocess
import longranger.cnv.contig_manager as contig_manager

__MRO__ = """
stage CALCULATE_METRICS(
    in bam[] merged_bams,
    in string[] node_ids,
    in string reference_path,
    out json metrics,
    src py "stages/calculate_metrics",
) split using(
    in bam merged_bam,
    in string node_id,
    out bw results,
)
"""

def split(args):
    return {'chunks': [{'subset_bam': bam, 'node_id': n, '__mem_gb': 6 * 16, '__threads': 16}
                       for bam, n in zip(args.merged_bams, args.node_ids)],
            'join': {'__mem_gb': 16}}


def main(args, outs):
    """
    Produces a bigWig format depth vector
    """
    cmd = ['bamCoverage', '-b', args.subset_bam, '-o', outs.results,
           '-bs', '1', '--ignoreDuplicates', '-p', str(args.__threads)]
    subprocess.check_call(cmd)


def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    ref = contig_manager.contig_manager(args.reference_path)
    # I am doing this in a hacky wiggletools way because it is fast
    results = collections.defaultdict(dict)
    for c, d in zip(chunk_outs, chunk_defs):
        # whole genome metrics
        for metric in ['meanI', 'varI', 'CVI']:
            r, _ = subprocess.Popen(['wiggletools', metric, c.results], stdout=subprocess.PIPE).communicate()
            results[d.node_id][metric] = float(r.rstrip())
        # per chromosome
        per_chrom = collections.defaultdict(dict)
        tot_breadth = 0
        for chrom, chrom_length in ref.contig_lengths.iteritems():
            for metric in ['meanI', 'varI', 'CVI']:
                cmd = 'wiggletools seek {} 0 {} {} | wiggletools {} - '.format(chrom, chrom_length,
                                                                               c.results, metric)
                r, _ = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).communicate()
                per_chrom[metric][chrom] = float(r.rstrip())
            # calculate breadth, parsing in python
            cmd = 'wiggletools seek {} 0 {} {} | wiggletools unit - '.format(chrom, chrom_length, c.results)
            r, _ = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).communicate()
            r = r.split('\n')[:-1]
            breadth = 0
            for l in r:
                try:
                    chrom, start, stop, _ = l.split()
                except:
                    continue
                breadth += int(stop) - int(start)
            per_chrom['breadth'][chrom] = breadth
            per_chrom['breadth_frac'][chrom] = 1.0 * breadth / chrom_length
            tot_breadth += breadth
        results[d.node_id]['per_chromosome'] = per_chrom
        results[d.node_id]['genome_breadth'] = tot_breadth
        results[d.node_id]['breadth_fraction'] = 1.0 * tot_breadth / sum(ref.contig_lengths.values())
    with open(outs.metrics, 'w') as outf:
        json.dump(results, outf)
