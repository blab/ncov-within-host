'''
Takes in a tsv file of metadata of household pairs.
Chooses "random" pairs with similar distribution of Ct, # of variants, and frequency of variants.
Combines household pairs & random pairs metadata one df.
'''
import argparse
import json
import pandas as pd
import numpy as np
from scipy import stats

def load_metadata(file):
    '''
    Loads metadata tsv as df.
    '''
    with open(file) as tfile:
        df = pd.read_csv(tfile, sep = '\t')
        df['nwgc_id'] = df.nwgc_id.astype('str')
        if 'origin' in df.columns:
            df.loc[df.origin=='sfs', 'origin'] = 'scan'
    return df

def load_snvs(file):
    '''
    Loads SNVs dictionary
    '''
    with open(file) as jfile:
        snvs = json.load(jfile)
    return snvs

def label_pairs(df):
    '''
    Adds pair type & numbers pairs in pairs.
    '''
    n_pairs = int(len(df)/2)
    df['pair'] = [pair for num in range(1, n_pairs + 1) for pair in [num]*2]
    df['pair_type'] = 'household'
    return df

def add_variants(df, snvs):
    '''
    Adds total number of SNV's to dataframe.
    '''
    n_snvs = []
    for sample in df['nwgc_id']:
        n_snvs.append(snvs[sample]['total'])
    df['n_snvs'] = n_snvs
    return df

def create_snvs_df(metadata, snvs):
    '''
    Constructs df for analyses
    '''
    positions = []
    variants = []
    frequencies = []
    ids = []
    for sample in metadata['nwgc_id']:
        positions.extend(snvs[sample]['position'])
        variants.extend(snvs[sample]['variant'])
        frequencies.extend(snvs[sample]['frequency'])
        ids.extend([sample]*len(snvs[sample]['position']))
    df = pd.DataFrame()
    df['nwgc_id'] = ids
    df['position'] = positions
    df['variant'] = variants
    df['frequency'] = frequencies
    df = df.merge(metadata, on='nwgc_id', how='left')
    return df

def choose_random_pairs(metadata, pairs, origin):
    '''
    Chooses random pairs as controls for household pairs.
    '''
    max_ct = max(pairs['avg_ct'])
    min_ct = min(pairs['avg_ct'])
    ids = list(pairs['nwgc_id'])
    n_pairs = int(len(pairs)/2)

    filtered = metadata[(metadata.avg_ct <= max_ct) & (metadata.avg_ct >= min_ct) & (metadata.origin == origin) & (~metadata.nwgc_id.isin(ids))]
    chosen = np.random.choice(filtered.nwgc_id, size = (7, 2), replace=False)
    chosen_list = chosen.flatten().tolist()
    df = metadata[metadata.nwgc_id.isin(chosen_list)]
    df = df.dropna(axis=1, how='all')

    df['pair'] = [pair for num in range(n_pairs + 1,n_pairs*2 + 1,1) for pair in [num]*2 ]
    df['pair_type'] = 'random'
    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Makes JSON of iSNVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True, help='metadata tsv')
    parser.add_argument('--pairs', type=str, required=True, help = 'tsv with household pairs')
    parser.add_argument('--snvs', type=str, required=True, help = 'json of snvs')
    parser.add_argument('--origin', type=str, required=True, help = 'sample origin')
    parser.add_argument('--output', type=str, required=True, help = 'location of output tsv')
    args = parser.parse_args()

    # Loads datasets
    metadata = load_metadata(args.metadata)
    pairs = load_metadata(args.pairs)
    snvs = load_snvs(args.snvs)

    # Labels household pairs & adds # of variants
    pairs = label_pairs(pairs)
    pairs = add_variants(pairs, snvs)

    # Creates snvs df for household pairs
    pairs_snvs = create_snvs_df(pairs, snvs)

    # Chooses random pairs with same distribution of Ct, variants, and frequency of variants as household pairs
    finished = False
    while not finished:
        rand_pairs = choose_random_pairs(metadata, pairs, args.origin)
        if stats.ranksums(rand_pairs['avg_ct'], pairs['avg_ct']).pvalue > 0.05:
            rand_pairs = add_variants(rand_pairs, snvs)
            if stats.ranksums(pairs['n_snvs'], rand_pairs['n_snvs']).pvalue > 0.05:
                rand_pairs_snvs = create_snvs_df(rand_pairs, snvs)
                if stats.ranksums(pairs_snvs['frequency'], rand_pairs_snvs['frequency']).pvalue > 0.05:
                    finished = True

    # Adds random pairs to household pairs dataframe
    all_pairs = pd.concat([pairs, rand_pairs],join='inner', ignore_index=True)

    # Saves final df
    with open(args.output, 'w') as f:
        all_pairs.to_csv(f, sep = '\t', index=False)
