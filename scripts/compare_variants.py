'''
Motivation: Needed a way to compare variants called on the same sample with
different pipelines.

Action: Compares variants from two different dataframes. Each dataframe must
contain the same samples. Output is a pdf  plotting frequency of variant in each
dataset + tsv of variants, their frequency, and coverage for each sample.

Flags:
--variants, dataframes containing variants
--output, path to output folder
--no-indel, include to only compare iSNVs and not indels.
--label, the names of each pipeline. Used in plots and tsv column headers.
'''

import argparse
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def load_df(file):
    '''
    Loads variants tsv as df.
    Creates snvs column
    '''
    with open(file) as tfile:
        df = pd.read_csv(tfile, sep = '\t')
        if 'nwgc_id' in df.columns:
            df['nwgc_id'] = df.nwgc_id.astype('str')
        if 'position' and 'variant' in df.columns:
            df['snv'] = df.loc[:,['reference', 'position', 'variant']].apply(
                lambda x: ''.join(x.dropna().astype(str)),
                axis=1
            )
    return df

def filter_frequencies(df, cutoff):
    '''
    Filters variants to above a cutoff frequency
    '''
    return df[df.frequency > cutoff].reset_index(drop=True)

def drop_indel(df):
    '''
    Drops indel variants
    '''
    return df[df.variant.str.len() == 1].reset_index(drop=True)

def compare_frequencies(df1, df2):
    '''
    Creates dataframe containing the frequency of each sample variant in each pipeline
    '''
    ids = []
    freq1 = []
    freq2 = []
    variants = []
    coverage1 = []
    coverage2 = []
    for nwgc_id in df1['nwgc_id'].unique():
        snvs1 = df1.loc[df1.nwgc_id == nwgc_id, 'snv']
        snvs2 = df2.loc[df2.nwgc_id == nwgc_id, 'snv']
        snvs = snvs1.append(snvs2).unique()
        for snv in snvs:
            if snv in snvs1.values:
                freq1.append(df1.loc[(df1.snv==snv) & (df1.nwgc_id==nwgc_id), 'frequency'].values[0])
                coverage1.append(df1.loc[(df1.snv==snv) & (df1.nwgc_id==nwgc_id), 'coverage'].values[0])
            else:
                freq1.append(0)
                coverage1.append(0)
            if snv in snvs2.values:
                freq2.append(df2.loc[(df2.snv==snv) & (df2.nwgc_id==nwgc_id), 'frequency'].values[0])
                coverage2.append(df2.loc[(df2.snv==snv) & (df2.nwgc_id==nwgc_id), 'coverage'].values[0])
            else:
                freq2.append(0)
                coverage2.append(0)
            ids.append(nwgc_id)
            variants.append(snv)

    frequencies = pd.DataFrame({ 'nwgc_id':ids, 'freq1':freq1, 'freq2':freq2, 'snv':variants, 'coverage1':coverage1, 'coverage2':coverage2})
    return frequencies

def plot_frequencies(df, label, output):
    '''
    Plots variant frequency in pipeline 1 vs. pipeline 2
    '''
    fig, (ax1, ax2) = plt.subplots(1,2)

    sns.scatterplot(data = df, x='freq1', y='freq2', ax=ax1)
    sns.scatterplot(data = df, x='freq1', y='freq2', ax=ax2)

    ax1.set_ylabel("Freq in " + label[1])
    ax2.set_ylabel("Freq in " + label[1])
    ax1.set_xlabel("Freq in " + label[0])
    ax2.set_xlabel("Freq in " + label[0])
    ax1.set_aspect('equal', adjustable='box')
    ax1.set_xlim(-0.05, 1)
    ax1.set_ylim(-0.05, 1)
    ax1.set_xticks(np.arange(0, 1.1, 0.1))
    ax1.set_yticks(np.arange(0, 1.1, 0.1))
    ax2.set_aspect('equal', adjustable='box')
    ax2.set_xlim(-0.005, 0.1)
    ax2.set_ylim(-0.005, 0.1)
    ax2.set_xticks(np.arange(0, 0.12, 0.02))
    ax2.set_yticks(np.arange(0, 0.12, 0.02))

    plt.tight_layout()
    return plt.savefig(fname=output, dpi=300)

def label_columns(df, label):
    '''
    Labels tsv columns with appropriate pipelines
    '''
    name_dict = {'freq1' : 'freq_' + label[0],
                'freq2' : 'freq_' + label[1],
                'coverage1' : 'coverage_' + label[0],
                'coverage2' : 'coverage_' + label[1]}
    return df.rename(columns=name_dict)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Makes JSON of iSNVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--variants', nargs = '+', type=str, required=True, help='path to input files')
    parser.add_argument('--cutoff', type=float, default=0.01, help='variant frequency cutoff')
    parser.add_argument('--no-indel', action='store_true', help = 'flag to only compare snvs, not indels')
    parser.add_argument('--label', nargs='+', type=str, required=True, help = 'pipeline labels')
    parser.add_argument('--output', type=str, required=True, help = 'location of output folder')
    args = parser.parse_args()

    # Loads SNVs dataframes
    snvs_A = load_df(args.variants[0])
    snvs_B = load_df(args.variants[1])

    # Filters to frequency
    filtered_A = filter_frequencies(snvs_A, args.cutoff)
    filtered_B = filter_frequencies(snvs_B, args.cutoff)

    # Removes indels
    if args.no_indel == True:
        filtered_A = drop_indel(filtered_A)
        filtered_B = drop_indel(filtered_B)

    # Creates df with frequency & coverge of each variant in each dataset
    compared = compare_frequencies(filtered_A, filtered_B)

    # Make ouput directories
    plots_path = args.output + 'plots/'

    if not os.path.isdir(args.output):
        os.mkdir(args.output)
    if not os.path.isdir(plots_path):
        os.mkdir(plots_path)

    # Plots variants
    for nwgc_id in compared.nwgc_id.unique():
        to_plot = compared[compared.nwgc_id == nwgc_id]
        plot_path = plots_path + str(nwgc_id) + '_freqs.pdf'
        plot_frequencies(to_plot, args.label, plot_path)

    # Labels columns in DataFrame
    labelled = label_columns(compared, args.label)

    # Write out tsv
    tsv_path = args.output + 'variant_frequencies.tsv'
    with open(tsv_path, 'w') as f:
        labelled.to_csv(f, sep = '\t', index=False)
