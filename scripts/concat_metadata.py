'''
Combines metadata from GISAID, Metabase, WA-DoH, and the strain-to-nwgc-id file.
'''

import argparse
import pandas as pd
from Bio import SeqIO
import numpy as np

def load_metadata(file, sep):
    '''
    Loads metadata tsv as df.
    '''
    with open(file) as tfile:
        metadata = pd.read_csv(tfile, sep = sep)
        if 'nwgc_id' in metadata.columns:
            metadata['nwgc_id'] = metadata.nwgc_id.astype('str')
    return metadata

def prep_global(df):
    '''
    Selects only necessary columns from global metadata.
    '''
    df1 = df[df['submitting_lab'].str.contains('Seattle Flu Study', regex=False)]
    df2 = df1[['strain', 'date', 'location', 'country_exposure', 'region_exposure', 'division_exposure']]
    return df2

def merge(df1, df2, on):
    '''
    Left merges df.
    '''
    df = df1.merge(df2, how='left', on=on)
    return df

def coalesce_address_symptom_ct(df):
    '''
    Left merges df + coalesces columns where the values are the same.
    '''
    df['address_identifier_x'].update(df['address_identifier_y'])
    df['symptom_onset_x'].update(df['symptom_onset_y'])
    df['avg_ct_x'].update(df['avg_ct_y'])
    df.rename(columns={"address_identifier_x": "address_identifier", "symptom_onset_x" : "symptom_onset", "avg_ct_x" : "avg_ct"}, inplace=True)
    df.drop(['address_identifier_y', 'symptom_onset_y', "avg_ct_y"], axis=1, inplace=True)
    return df


def add_reference(df):
    '''
    Adds column containing reference strain for all samples.
    '''
    df['reference'] = df['strain']
    for index, row in df[df.reference.isna()].iterrows():
        if row['nwgc_id'] in df[df.index != index]['nwgc_id'].unique():
            ref = df[(df.index != index) & (df.nwgc_id == row['nwgc_id'])]['strain'].dropna()
        elif row['sample_barcode'] in df[df.index != index]['sample_barcode'].unique():
            ref = df.loc[(df.index != index) & (df.sample_barcode == row['sample_barcode'])]['strain'].dropna()
        elif row['phl_accession'] in df[df.index != index]['phl_accession'].unique():
            ref = df.loc[(df.index != index) & (df.phl_accession == row['phl_accession'])]['strain'].dropna()
        else:
            ref = pd.Series(np.nan)
        df.loc[index, 'reference'] = ref.item()
    return df

def filter_metadata(fasta, df):
    '''
    Removes samples from metadata not in fasta file.
    '''
    genomes = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    strains = [strain for strain in list(df['strain']) if strain not in genomes]
    filtered = df[~df.strain.isin(strains)]
    return filtered

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Makes JSON of iSNVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--strains', type=str, required=True, help = 'strain to nwgc_id csv')
    parser.add_argument('--wadoh', type=str, required=True, help = 'WA-DoH metadata')
    parser.add_argument('--global-metadata', type=str, required=True, help = 'global metadata')
    parser.add_argument('--metabase', type=str, required=True, help = 'Metabase metadata')
    parser.add_argument('--wdrs', type=str, required=True, help = 'WDRS metadata')
    parser.add_argument('--genomes', type=str, required=True, help = 'Fasta')
    parser.add_argument('--output', type=str, required=True, help = 'location of output tsv')
    args = parser.parse_args()

    # Load all metadata to merge
    strains = load_metadata(args.strains, '\t')
    wadoh = load_metadata(args.wadoh, '\t')
    global_metadata = load_metadata(args.global_metadata, '\t')
    wdrs = load_metadata(args.wdrs, '\t')
    metabase = load_metadata(args.metabase, '\t')

    # Prep metadata for merging
    global_prepped = prep_global(global_metadata)

    # Merge dataframes
    dfA = merge(strains, wadoh, 'nwgc_id')
    dfB = merge(dfA, global_prepped, 'strain')
    dfC = merge(dfB, metabase, 'sample_barcode')
    dfD = merge(dfC, wdrs, 'strain')
    dfE = coalesce_address_symptom_ct(dfD)
    dfF = add_reference(dfE)

    # Save final df
    with open(args.output, 'w') as f:
        dfF.to_csv(f, sep = '\t', index=False)
