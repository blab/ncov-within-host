'''
This script combines WA-DoH metadata into a single file.
'''

import argparse
import pandas as pd

def merge(flist):
    '''
    Merges metadata.
    '''
    counter = 0
    for file in flist:
        with open(file) as tfile:
            metadata = pd.read_csv(tfile, sep = '\t')
            if counter == 0:
                df = metadata
                counter += 1
            else:
                df = df.merge(metadata, how = 'outer')
    return df

def clean(file):
    '''
    Cleans metadata
    '''
    df1 = file.drop(columns = ['sex', 'investigator_id', 'collection_date', 'county'])
    df2 = df1.dropna(subset=['nwgc_id'])
    df3 = df2.astype({'nwgc_id': 'int32'})
    df4 = df3.dropna(thresh=2)
    df5 = df4.loc[df4.primers.isna(), :]
    df6 = df4.loc[df4.primers.notna(), :]
    df5.set_index('nwgc_id', inplace=True)
    df6.set_index('nwgc_id', inplace=True)
    df7 = df5.combine_first(df6).reset_index()
    return df7

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Makes JSON of iSNVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, nargs='+', required=True, help = 'list of tsv files containing metadata')
    parser.add_argument('--output', type=str, required=True, help = 'location of output tsv')
    args = parser.parse_args()

    # Merge metadata
    df = merge(args.metadata)

    # Clean metadata
    df_clean = clean(df)

    # Write out metadata
    with open(args.output, 'w') as f:
        df_clean.to_csv(f, sep = '\t', index=False)
