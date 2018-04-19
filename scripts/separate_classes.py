import argparse
import os
import pandas as pd


if __name__ == '__main__':
    args = argparse.ArgumentParser()
    args.add_argument('classification_file')
    args.add_argument('--outdir', default=os.getcwd())
    args = args.parse_args()
    df = pd.read_csv(args.classification_file)
    df['fixed_sample'] = [x.split('-')[1] for x in df['sample']]
    for (cell_type, sample), d in df.groupby(['cell_type', 'fixed_sample']):
        if cell_type == 'noisy':
            continue
        with open('{}-{}.txt'.format(sample, cell_type), 'w') as outf:
            for x in d.barcode:
                outf.write(x + '\n')