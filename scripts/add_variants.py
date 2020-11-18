import argparse
import json
import pandas as pd

def load_metadata(file):
    '''
    Loads metadata tsv as df.
    '''
    with open(file) as tfile:
        df = pd.read_csv(tfile, sep = '\t')
        df['nwgc_id'] = df.nwgc_id.astype('str')
        if 'origin' in df.columns:
            df.loc[df.origin=='sfs', 'origin'] = 'scan'
        df = df.sort_values(by=['address_identifier'], axis=0, ignore_index=True)
    return df

def load_snvs(file):
    '''
    Loads SNVs dictionary
    '''
    with open(file) as jfile:
        snvs = json.load(jfile)
    return snvs

def add_variants(df, snvs):
    '''
    Adds total number of SNV's to dataframe.
    '''
    n_snvs = []
    for sample in df['nwgc_id']:
        n_snvs.append(snvs[sample]['total'])
    df['n_snvs'] = n_snvs
    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Makes JSON of iSNVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True, help='metadata tsv')
    parser.add_argument('--snvs', type=str, required=True, help = 'json of snvs')
    parser.add_argument('--output', type=str, required=True, help = 'location of output tsv')
    args = parser.parse_args()

    # Loads datasets
    metadata = load_metadata(args.metadata)
    snvs = load_snvs(args.snvs)

    # Adds # of variants & checks that all samples have >0 variants.
    w_nsnvs = add_variants(metadata, snvs)

    # Saves final df
    with open(args.output, 'w') as f:
        w_nsnvs.to_csv(f, sep = '\t', index=False)
