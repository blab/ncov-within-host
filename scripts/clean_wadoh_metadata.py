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

def clean(df):
    '''
    Cleans metadata
    '''
    df = df.drop(columns = ['sex', 'collection_date', 'county', 'investigator_id'])
    df = df.dropna(subset=['nwgc_id'])
    df = df.astype({'nwgc_id': 'int32'})
    df = df.reset_index(drop=True)
    return df

def dedup(df):
    df['duplicate'] = df.duplicated(subset=['nwgc_id'], keep=False)
    for index, row in df.iterrows():
        if row.duplicate == True:
            duplicate = df.loc[(df.nwgc_id == row.nwgc_id) & (df.index != index), :]
            df.iloc[index, :] = row.fillna(duplicate.iloc[0])
    df.drop(columns=['duplicate'],axis=1, inplace=True)
    df.drop_duplicates(inplace=True, ignore_index = True)
    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Makes JSON of iSNVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, nargs='+', required=True, help = 'list of tsv files containing metadata')
    parser.add_argument('--output', type=str, required=True, help = 'location of output tsv')
    args = parser.parse_args()

    # Merge metadata
    merged = merge(args.metadata)

    # Clean metadata
    cleaned = clean(merged)

    # Coalesce duplicate rows
    deduped = dedup(cleaned)

    # Write out metadata
    with open(args.output, 'w') as f:
        deduped.to_csv(f, sep = '\t', index=False)
