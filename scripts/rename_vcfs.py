import argparse
import pandas as pd
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Makes JSON of iSNVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True, help='tsv with nwgc_id and id')
    parser.add_argument('--directory', type=str, required=True, help = 'folder with files')
    parser.add_argument('--ext', type=str, required=True, help = 'extension of files to rename')
    parser.add_argument('--batch', type=str, required=True, help = 'batch name')
    args = parser.parse_args()

    with open(args.metadata) as tfile:
        metadata = pd.read_csv(tfile, sep = '\t')
        if 'nwgc_id' in metadata.columns:
            metadata['nwgc_id'] = metadata.nwgc_id.astype('str')

    for path in os.listdir(args.directory):
        file  = path.split('/')[-1]
        nwgc_id = file.split('.')[0]
        if nwgc_id in metadata['nwgc_id'].unique():
            os.rename(args.directory + path, args.directory + args.batch + '_' + nwgc_id + '.' + args.ext)
    #        if len(metadata.loc[metadata.nwgc_id == nwgc_id]['id']) == 1:
    #            id = metadata.loc[metadata.nwgc_id == nwgc_id]['id'].item()
    #            os.rename(args.directory + nwgc_id + '.' + args.ext, args.directory + id + '.' + args.ext)
    #        else:
    #            print(nwgc_id  + ' has multiples.')

        #else:
        #    print(nwgc_id + ' not in metadata.')
