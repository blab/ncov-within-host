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
    reads = []

    vcf = allel.read_vcf(file, fields=['variants/POS',
                                       'variants/REF',
                                       'variants/ALT',
                                       'calldata/DP',
                                       'calldata/AD'],
                                       types={'calldata/DP': 'i4', 'calldata/AD':'i4'})
    if vcf is not None:
        for i in range(len(vcf['variants/POS'])):
            position.append(int(vcf['variants/POS'][i]))
            variant.append(vcf['variants/ALT'][i,0])
            reference.append(vcf['variants/REF'][i])
            coverage.append(int(vcf['calldata/DP'][i,0]))
            reads.append(int(vcf['calldata/AD'][i,0,0]))
            freq = float(vcf['calldata/AD'][i,0,0]/vcf['calldata/DP'][i,0])
            frequency.append(freq)

    df = pd.DataFrame({ 'position':position,
                        'variant':variant,
                        'reference':reference,
                        'frequency':frequency,
                        'coverage':coverage,
                        'reads':reads})
    return df

def multiplex_mapping(vcfs, metadata):
    '''
    Constructs df containing all variants from a list of vcfs
    '''
    df = pd.DataFrame()
    for file in vcfs:
        path = file.split('/')
        id = path[-1].split('.')[0]
        if metadata['id'].str.contains(id).any():
            sample = map_variants(file)
            sample['id'] = id
            df = df.append(sample)

    df = df.merge(metadata[['strain', 'id', 'avg_ct', 'nwgc_id']], how = 'left')
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
    parser.add_argument('--vcf', nargs = '+', type=str, required=True, help = 'list of vcf files')
    parser.add_argument('--output', type=str, required=True, help = 'location of output tsv')
    args = parser.parse_args()

    # Load metadata
    metadata = load_df(args.metadata)

    # Map variants to df
    variants = multiplex_mapping(args.vcf, metadata)

    # Reorder columns
    columns_order = ['strain', 'id', 'position', 'variant', 'reference']
    ordered = reorder_cols(variants, columns_order)

    # Save variants
    with open(args.output, 'w') as f:
        ordered.to_csv(f, sep = '\t', index=False)
