'''
Combines metadata from GISAID, Metabase, WA-DoH, and the strain-to-nwgc-id file.
'''

import argparse
import pandas as pd
from Bio import SeqIO

def load_metadata(file, sep):
    '''
    Loads metadata tsv as df.
    '''
    with open(file) as tfile:
        metadata = pd.read_csv(tfile, sep = sep)
    return metadata

def prep_strains(df):
    '''
    Eliminates extraneous text from strain name.
    '''
    df['strain'] = df['strain'].str.replace('^.*?(/)', '', regex=True)
    return df

def prep_global(df):
    '''
    Selects only necessary columns from global metadata.
    '''
    df1 = df[df['submitting_lab'] == 'Seattle Flu Study']
    df2 = df1[['strain', 'date', 'location', 'country_exposure', 'region_exposure', 'division_exposure']]
    return df2

def merge(df1, df2, on):
    '''
    Left merges df.
    '''
    df = df1.merge(df2, how='left', on=on)
    return df

def merge_metabase(df, metadata):
    '''
    Merges metabase df with other df.
    '''
    df['individual_identifier'] = ''
    df['clia_barcode'] = ''
    df['age'] = ''
    df['puma'] = ''
    df['address_identifier'] = ''
    df['symptom_onset'] = ''
    for i in df.index:
        for j in metadata.index:
            if str(df.nwgc_id[i]) in metadata.nwgc_id[j]:
                df.loc[i, ['individual_identifier', 'clia_barcode', 'age', 'puma', 'address_identifier', 'symptom_onset', 'avg_ct']] = metadata.loc[j, ['individual_identifier', 'clia_barcode', 'age', 'puma', 'address_identifier', 'symptom_onset', 'avg_ct']]
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
    parser.add_argument('--genomes', type=str, required=True, help = 'Fasta')
    parser.add_argument('--output', type=str, required=True, help = 'location of output tsv')
    args = parser.parse_args()

    # Load all metadata to merge
    strains = load_metadata(args.strains, ',')
    wadoh = load_metadata(args.wadoh, '\t')
    global_metadata = load_metadata(args.global_metadata, '\t')
    metabase = load_metadata(args.metabase, '\t')

    # Prep metadata for merging
    strains_prepped = prep_strains(strains)
    global_prepped = prep_global(global_metadata)

    # Merge dataframes
    dfA = merge(strains_prepped, wadoh, 'nwgc_id')
    dfB = merge(dfA, global_prepped, 'strain')
    dfC = merge_metabase(dfB, metabase)
    dfD = filter_metadata(args.genomes, dfC)

    # Save final df
    with open(args.output, 'w') as f:
        dfD.to_csv(f, sep = '\t', index=False)
