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
    return metadata

def df_to_plot(snvs, metadata):
    '''
    Constructs dataframe to plot.
    '''
    samples = [key for key in snvs.keys()]
    freqs = [frequency for frequency in snvs[samples[0]]['total'].keys()]
    df = pd.DataFrame(np.repeat(samples, len(freqs), axis=0), columns=['nwgc_id'])
    df['nwgc_id'] = df['nwgc_id'].astype('i4')

    cycle_freqs = cycle(freqs)
    df['freq'] = [next(cycle_freqs) for i in range(len(df))]

    df = df.merge(metadata[['nwgc_id', 'avg_ct', 'origin']], how = 'left')

    for sample in samples:
        for freq in freqs:
            df.loc[(df.nwgc_id == int(sample)) & (df.freq == freq), 'isnvs'] = snvs[sample]['total'][freq]
    df.dropna(subset = ['avg_ct'], inplace=True)
    return df

def plot(df, output):
    '''
    Plots iSNVs vs. Ct at minor allele frequencies specied by maf.
    '''
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(8, 8), facecolor='white')

    for freq in sorted(set(df['freq'])):
        x_s = df.loc[(df.freq == freq) & (df.origin == 'sfs'), 'avg_ct']
        y_s = df.loc[(df.freq == freq) & (df.origin == 'sfs'), 'isnvs']
        ax1.scatter(x_s, y_s, label=freq, edgecolors='none')
        ax2.scatter(x_s, y_s, label='', edgecolors='none')

        x_w = df.loc[(df.freq == freq) & (df.origin == 'wadoh'), 'avg_ct']
        y_w = df.loc[(df.freq == freq) & (df.origin == 'wadoh'), 'isnvs']
        ax3.scatter(x_w, y_w, label='', edgecolors='none')
        ax4.scatter(x_w, y_w, label='', edgecolors='none')

    ax1.set_title('SCAN')
    ax2.set_title('SCAN: Zoomed in')
    ax3.set_title('WA-DoH')
    ax4.set_title('WA-DOH: Zoomed in')

    ax2.set_ylim(bottom=-1, top = 41)
    ax2.set_xlim(left=15.5, right = 25.5)
    ax4.set_ylim(bottom=-1, top = 51)
    ax4.set_xlim(left=15, right = 25.5)

    ax1.set_xlabel('Ct')
    ax2.set_xlabel('Ct')
    ax3.set_xlabel('Ct')
    ax4.set_xlabel('Ct')

    ax1.set_ylabel('iSNVs per sample')
    ax2.set_ylabel('iSNVs per sample')
    ax3.set_ylabel('iSNVs per sample')
    ax4.set_ylabel('iSNVs per sample')

    fig.legend(title = 'Minor variant\nfrequencies', loc = 'center', bbox_to_anchor=(1.05, 0.5),
               labels = {'0.01':'mvf=0.01', '0.02':'mvf=0.02', '0.05':'mvf=0.05'})
    fig.tight_layout()
    return plt.savefig(output, dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Makes JSON of iSNVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--snvs', type=str, required=True, help='iSNVs json')
    parser.add_argument('--metadata', type=str, required=True, help = 'tsv file containing metadata')
    parser.add_argument('--output', type=str, required=True, help = 'location of output json')
    args = parser.parse_args()

    # Load iSNVs
    snvs = load_snvs(args.snvs)

    # Load metadata
    metadata = load_metadata(args.metadata)

    # Construct df to plots
    df = df_to_plot(snvs,  metadata)

    # Save plot to file
    plot(df, args.output)
