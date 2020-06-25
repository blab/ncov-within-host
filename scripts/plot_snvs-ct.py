'''
This script plots Ct vs. number of iSNVs per SARS-CoV-2 genomes at various minor allele frequencies.
'''
import argparse
import json
import pandas as pd
import numpy as np
from itertools import cycle
import matplotlib.pyplot as plt

def load_snvs(file):
    '''
    Loads JSON of SNVs as dictionary.
    '''
    with open(file) as jfile:
        snvs = json.load(jfile)
    return snvs

def load_metadata(file):
    '''
    Loads metadata tsv as df.
    '''
    with open(file) as tfile:
        metadata = pd.read_csv(tfile, sep = '\t')
        metadata.rename(columns={'avg_hcov19_crt':'ct'}, inplace=True)
        sample_ids = metadata['nwgc_id'].apply(lambda x:
                                      np.fromstring(
                                          x.replace('[', '')
                                          .replace(']', ''), sep =',', dtype=int))
        metadata = metadata.assign(nwgc_id = sample_ids)
    return metadata

def df_to_plot(snvs, maf, metadata):
    '''
    Constructs dataframe to plot.
    '''
    samples = [key for key in snvs.keys()]
    df = pd.DataFrame(np.repeat(samples, len(maf), axis=0), columns=['nwgc_id'])

    freqs = cycle(maf)
    df['freq'] = [next(freqs) for freq in range(len(df))]

    for sample in samples:
        for i in metadata.nwgc_id.index:
            if int(sample) in metadata.nwgc_id[i]:
                ct = metadata.iloc[i, 6]
                df.loc[df.nwgc_id == sample, 'ct'] = ct
        for freq in snvs[sample].keys():
            df.loc[(df.nwgc_id == sample) & (df.freq == float(freq)), 'isnvs'] = snvs[sample][freq]['total']
    df.dropna(inplace=True)
    return df

def plot(maf, df, output):
    '''
    Plots iSNVs vs. Ct at minor allele frequencies specied by maf.
    '''
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(8, 4), facecolor='white')
    for freq in maf:
        x = df.loc[df.freq == freq, 'ct']
        y = df.loc[df.freq == freq, 'isnvs']
        ax1.scatter(x, y, label=freq, edgecolors='none')
        ax2.scatter(x, y, label='', edgecolors='none')
    ax2.set_title('Zoomed in')
    ax2.set_ylim(bottom=-1, top = 41)
    ax2.set_xlim(left=15.5, right = 25.5)
    ax1.set_xlabel('Ct')
    ax2.set_xlabel('Ct')
    ax1.set_ylabel('iSNVs per sample')
    ax2.set_ylabel('iSNVs per sample')
    fig.legend(title = 'Minor allele\nfrequencies', loc = 'center', bbox_to_anchor=(1.05, 0.5))
    fig.tight_layout()
    return plt.savefig(output, dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Makes JSON of iSNVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--snvs', type=str, required=True, help='iSNVs json')
    parser.add_argument('--metadata', type=str, required=True, help = 'tsv file containing metadata')
    parser.add_argument('--maf', nargs='+', type=float, required=True, help = 'list of minor allele frequencies')
    parser.add_argument('--output', type=str, required=True, help = 'location of output json')
    args = parser.parse_args()

    # Load iSNVs
    snvs = load_snvs(args.snvs)

    # Load metadata
    metadata = load_metadata(args.metadata)

    # Construct df to plots
    df = df_to_plot(snvs, args.maf, metadata)

    # Save plot to file
    plot(args.maf, df, args.output)
