'''
Generates TV plots comparing SNVs per replicate set and for all replicates.
Flags:
 --variants, path to df with snvs/indels
 --replicates, path to tsv mapping replicates
 --cutoff, (optional)variant frequency cutoff
 --output, path to output directory
'''

import argparse
import pandas as pd
import os
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def load_variants(file):
    '''
    Loads variants tsv as df.
    Creates snvs column
    '''
    with open(file) as tfile:
        df = pd.read_csv(tfile, sep = '\t')
        if 'nwgc_id' in df.columns:
            df['nwgc_id'] = df.nwgc_id.astype('str')

        if 'position' and 'variant' in df.columns:
            df['position'] = pd.to_numeric(df.position, downcast='integer')
            df['snv'] = df.loc[:,['reference', 'position', 'variant']].apply(
                lambda x: ''.join(x.dropna().astype(str)),
                axis=1
            )
    return df

def load_replicates(file):
    '''
    Loads replicates tsv as df.
    '''
    with open(file) as tfile:
        df = pd.read_csv(tfile, sep = '\t')
    return df

def filter_frequencies(df, cutoff):
    '''
    Filters variants to above a cutoff frequency
    '''
    return df[df.frequency > cutoff].reset_index(drop=True)

def compare_frequencies(vDF, rDF, pileups):
    '''
    Creates dataframe containing the frequency of each sample variant in each pipeline
    '''
    replicates = []
    id1 = []
    id2 = []
    freq1 = []
    freq2 = []
    variants = []
    coverage1 = []
    coverage2 = []
    avg_ct = []
    reads1 = []
    reads2 = []
    for s1, s2, num in zip(rDF['id1'], rDF['id2'], rDF['replicate_no']):
        ct = vDF.loc[vDF.id == s1, 'avg_ct'].unique()[0]
        with open(pileups + s1 + '.pileup') as tfile:
            p1 = pd.read_csv(tfile, sep = '\t', header=None, nrows=29903)
        with open(pileups + s2 + '.pileup') as tfile:
            p2 = pd.read_csv(tfile, sep = '\t', header=None, nrows=29903)
        snvs1 = vDF.loc[vDF.id == s1, 'snv']
        snvs2 = vDF.loc[vDF.id == s2, 'snv']
        snvs = snvs1.append(snvs2).unique()
        positions = [int(snv[1:-1]) for snv in snvs]
        for snv, pos in zip(snvs, positions):
            if snv in snvs1.values:
                freq1.append(vDF.loc[(vDF.snv==snv) & (vDF.id==s1), 'frequency'].values[0])
                reads1.append(vDF.loc[(vDF.id==s1) & (vDF.snv==snv), 'reads'].values[0])
            else:
                freq1.append(0)
                reads1.append(0)
            if snv in snvs2.values:
                freq2.append(vDF.loc[(vDF.snv==snv) & (vDF.id==s2), 'frequency'].values[0])
                reads2.append(vDF.loc[(vDF.id==s2) & (vDF.snv==snv), 'reads'].values[0])
            else:
                freq2.append(0)
                reads2.append(0)
            coverage1.append(p1.at[pos-1,3])
            coverage2.append(p2.iat[pos-1,3])
            id1.append(s1)
            id2.append(s2)
            avg_ct.append(ct)
            replicates.append(num)
            variants.append(snv)

    frequencies = pd.DataFrame({ 'id1':id1, 'id2':id2, 'replicate_no':replicates, 'freq1':freq1, 'freq2':freq2, 'snv':variants, 'coverage1':coverage1, 'coverage2':coverage2, 'reads1':reads1, 'reads2':reads2, 'avg_ct':avg_ct})
    return frequencies

def plot_frequencies(df, output, type):
    '''
    Plots variant frequency in pipeline 1 vs. pipeline 2
    '''
    fig, (ax1, ax2) = plt.subplots(1,2)

    sns.scatterplot(data = df[(df.coverage1 >= 100) & (df.coverage2 >= 100)], x='freq1', y='freq2', ax=ax1)
    sns.scatterplot(data = df[(df.coverage1 >= 100) & (df.coverage2 >= 100)], x='freq1', y='freq2', ax=ax2)

    if type == 'sample':
        ax1.set_ylabel("Freq in " + df['id2'].unique())
        ax2.set_ylabel("Freq in " + df['id2'].unique())
        ax1.set_xlabel("Freq in " + df['id1'].unique())
        ax2.set_xlabel("Freq in " + df['id1'].unique())
    ax1.set_aspect('equal', adjustable='box')
    ax1.set_xlim(-0.05, 1.05)
    ax1.set_ylim(-0.05, 1.05)
    ax1.set_xticks(np.arange(0, 1.1, 0.1))
    ax1.set_yticks(np.arange(0, 1.1, 0.1))
    ax2.set_aspect('equal', adjustable='box')
    ax2.set_xlim(-0.005, 0.1)
    ax2.set_ylim(-0.005, 0.1)
    ax2.set_xticks(np.arange(0, 0.12, 0.02))
    ax2.set_yticks(np.arange(0, 0.12, 0.02))
    #ax1.set_xticklabels(ax1.get_xticklabels(), rotation='vertical')

    plt.tight_layout()
    plt.savefig(fname=output, dpi=300)
    return plt.close()

#def plot_histograms(df, output):
    '''
    Plots density histogram of snvs by Ct, coverage, frequency, # of reads.
    Colors by intersection snv vs. not.
    '''
#    fig,
#    return plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Makes JSON of iSNVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--variants', type=str, required=True, help='path to snvs df')
    parser.add_argument('--replicates', type=str, required=True, help='path to replicate mapping tsv')
    parser.add_argument('--pileups', type=str, required=True, help='path to directory with pileup files')
    parser.add_argument('--cutoff', type=float, default=0.01, help='variant frequency cutoff')
    parser.add_argument('--output', type=str, required=True, help = 'location of output folder')
    args = parser.parse_args()

    # Loads dataframes
    variants = load_variants(args.variants)
    replicates = load_replicates(args.replicates)

    # Filters to frequency
    filtered = filter_frequencies(variants, args.cutoff)

    # Creates df with frequency & coverage of each variant in each dataset
    compared = compare_frequencies(filtered, replicates, args.pileups)

    # Make ouput directories
    plots_path = args.output + 'plots/'

    if not os.path.isdir(args.output):
        os.mkdir(args.output)
    if not os.path.isdir(plots_path):
        os.mkdir(plots_path)

    # Plots variants
    for num in compared.replicate_no.unique():
        to_plot = compared[compared.replicate_no == num]
        plot_path = plots_path + 'rep_' + str(num) + '_freqs.pdf'
        plot_frequencies(to_plot, plot_path, 'sample')
    summary_path = args.output + 'freq1_freq2_summary.pdf'
    plot_frequencies(compared, summary_path, 'nope')

    #hist_path = args.output + 'replicates_histograms.pdf'
    #plot_histograms(compared, hist_path)

    # Write out tsv
    tsv_path = args.output + 'variant_frequencies.tsv'
    with open(tsv_path, 'w') as f:
        compared.to_csv(f, sep = '\t', index=False)
