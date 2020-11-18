'''
Takes in list of vcf files; outputs a tsv containing variants in relation to
Wuhan-Hu-1.

Inputs are:
    --vcfs, list of vcf files
    --output, location of output tsv
'''

import argparse
import os
os.environ["NUMEXPR_MAX_THREADS"]="272" # Prevents NUM_EXPr error that occurs when loading sci-kit allel
import pandas as pd
import allel

def load_df(file):
    '''
    Loads metadata tsv as df.
    '''
    with open(file) as tfile:
        df = pd.read_csv(tfile, sep = '\t')
        if 'nwgc_id' in df.columns:
            df['nwgc_id'] = df.nwgc_id.astype('str')
    return df

def map_variants(file):
    '''
    Maps variants in one vcf to df
    '''
    position = []
    variant = []
    reference = []
    frequency = []
    coverage = []

    vcf = allel.read_vcf(file, fields=['variants/POS',
                                       'variants/REF',
                                       'variants/ALT',
                                       'calldata/DP',
                                       'calldata/RD',
                                       'calldata/AD'],
                                       types={'calldata/DP': 'i4', 'calldata/AD':'i4'})
    if vcf is not None:
        for i in range(len(vcf['variants/POS'])):
            position.append(int(vcf['variants/POS'][i]))
            variant.append(vcf['variants/ALT'][i,0])
            reference.append(vcf['variants/REF'][i])
            coverage.append(int(vcf['calldata/DP'][i,0]))
            freq = float(vcf['calldata/AD'][i,0,0]/vcf['calldata/DP'][i,0])
            frequency.append(freq)

    df = pd.DataFrame({ 'position':position,
                        'variant':variant,
                        'reference':reference,
                        'frequency':frequency,
                        'coverage':coverage})
    return df

def multiplex_mapping(directory, metadata):
    '''
    Constructs df containing all variants from a list of vcfs
    '''
    df = pd.DataFrame()
    vcfs = os.listdir(directory)

    for file in vcfs:
        nwgc_id = file.split('.')[0]
        if metadata['nwgc_id'].str.contains(nwgc_id).any():
            sample = map_variants(directory + file)
            sample['nwgc_id'] = nwgc_id
            df = df.append(sample)

    df = df.merge(metadata[['strain', 'nwgc_id', 'avg_ct']], how = 'left')
    return df

def reorder_cols(df, cols):
    '''
    Reorders columns in cols to be first.
    '''
    return df[cols+list(set(df.columns).difference(cols))]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Makes JSON of iSNVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True, help='tsv with sample names and nwgc_id')
    parser.add_argument('--vcf', type=str, required=True, help = 'directory containing vcf files')
    parser.add_argument('--output', type=str, required=True, help = 'location of output tsv')
    args = parser.parse_args()

    # Load metadata
    metadata = load_df(args.metadata)

    # Map variants to df
    variants = multiplex_mapping(args.vcf, metadata)

    # Reorder columns
    columns_order = ['strain', 'nwgc_id', 'position', 'variant', 'reference']
    ordered = reorder_cols(variants, columns_order)

    # Save variants
    with open(args.output, 'w') as f:
        ordered.to_csv(f, sep = '\t', index=False)
