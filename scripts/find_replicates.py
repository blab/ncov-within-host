'''
Generates TSV of all replicate pairs.
'''

import argparse
import pandas as pd
import itertools

def load_df(file):
    '''
    Loads metadata tsv as df.
    '''
    with open(file) as tfile:
        df = pd.read_csv(tfile, sep = '\t')
        if 'nwgc_id' in df.columns:
            df['nwgc_id'] = df.nwgc_id.astype('str')
    return df

def find_duplicates(a):
    '''
    Returns duplicate values in 1D array.
    '''
    seen = {}
    dupes = []

    for x in a:
        if x not in seen:
            seen[x] = 1
        else:
            if seen[x] == 1:
                dupes.append(x)
            seen[x] += 1
    return dupes

def map_replicates(columns, df, exclude):
    '''
    Generates tsv with mapping of id for each replicate pair.
    '''
    id1 = []
    id2 = []
    matching_c = []
    matching_i = []
    for c in columns:
        dups = find_duplicates(df[c])
        for d in dups:
            ids = df.loc[df[c] == d, 'id']
            combs = list(itertools.combinations(ids, 2))
            for (a, b) in combs:
                id1.append(a)
                id2.append(b)
                matching_c.append(c)
                matching_i.append(d)
    reps = pd.DataFrame()
    reps['id1'] = id1
    reps['id2'] = id2
    reps['matching_variable'] = matching_c
    reps['matching_identifier'] = matching_i
    reps = reps[~reps.matching_identifier.isin(exclude)]
    dedup = reps.drop_duplicates(subset = ['id1', 'id2'])
    dedup['replicate_no'] = range(1, len(dedup) + 1)
    return dedup

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Makes JSON of iSNVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True, help='tsv file containing metadata')
    parser.add_argument('--match-on', nargs = '+', type=str, required=True, help = 'List of column names to match on')
    parser.add_argument('--exclude', nargs='+', type=str, default='[]', help= 'List of identifiers to exclude')
    parser.add_argument('--output', type=str, required=True, help = 'location of output tsv')
    args = parser.parse_args()

    # Load metadata
    metadata = load_df(args.metadata)

    # Generate replicate pairs
    replicates = map_replicates(args.match_on, metadata, args.exclude)

    # Write out tsv
    with open(args.output, 'w') as f:
        replicates.to_csv(f, sep = '\t', index=False)
