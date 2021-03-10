'''
In VCF files, SNVs are called relative to reference genome.
This generates a json with iSNVs relative to sample genome.

Inputs are:
    --metadata, a tsv file mapping strain to nwgc_id
    --sequences, fasta with consensus genome for each sample
    --vcfs, list of vcf files
    --output, location of json storing iSNVs
'''

import argparse
import pandas as pd
from Bio import SeqIO
import numpy as np
import os
os.environ["NUMEXPR_MAX_THREADS"]="272" # Prevents NUM_EXPr error that occurs when loading sci-kit allel
import allel
import json

def load_genomes(fasta):
    '''
    Loads fasta of genomes as dictionary.
    '''
    with open(fasta, 'r') as f:
        genomes = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
    return genomes

#def check_genomes(genomes, mfile):
#    '''
#    Checks that reference genome are available.
#    '''
#    with open(mfile, 'r') as f:
#        metadata = pd.read_csv(f, sep = '\t')
#    refs = set(metadata['reference'])
#    intersect = set(genomes.keys()).intersection(refs)
#    absent = [r for r in refs if r not in intersect]
#    print(absent)

def check_variants(file, ref, genomes):
    '''
    Opens vcf file and checks if SNVs based off reference genome are SNVs relative to sample consensus genome.
    If so, adds SNVs to dictionary supplied by mapping.
    '''
    mapping = {}
    mapping['position'] = []
    mapping['variant'] = []
    mapping['coverage'] = []
    mapping['frequency'] = []

    vcf = allel.read_vcf(file, fields=['variants/POS',
                                       'variants/REF',
                                       'variants/ALT',
                                       'calldata/DP',
                                       'calldata/RD',
                                       'calldata/AD',
                                       'calldata/RDF',
                                       'calldata/RDR',
                                       'calldata/ADF',
                                       'calldata/ADR'],
                                       types={'calldata/DP': 'i4', 'calldata/AD':'i4'})

    if vcf is not None:
        for i in range(len(vcf['variants/POS'])):
            pos = int(vcf['variants/POS'][i])

            if genomes[ref].seq[(pos-1)].upper() == vcf['variants/REF'][i]: # Minus 1 corrects for position numbering beginning at 0 for SeqIO records.
                mapping['position'].append(pos)
                mapping['variant'].append(vcf['variants/ALT'][i,0])
                mapping['coverage'].append(int(vcf['calldata/DP'][i,0]))
                mapping['frequency'].append(float(vcf['calldata/AD'][i,0,0]/vcf['calldata/DP'][i,0]))

            elif genomes[ref].seq[(pos-1)].upper() == vcf['variants/ALT'][i,0]: # Minus 1 corrects for position numbering beginning at 0 for SeqIO records.
                if (vcf['calldata/RD'][i,0])/(vcf['calldata/DP'][i,0]) >= 0.01: # Checks that reference variant meets cutoff
                    mapping['position'].append(pos)
                    mapping['variant'].append(vcf['variants/REF'][i])
                    mapping['coverage'].append(int(vcf['calldata/DP'][i,0]))
                    mapping['frequency'].append(float(vcf['calldata/RD'][i,0]/vcf['calldata/DP'][i,0]))

            elif genomes[ref].seq[(pos-1)].upper() == 'N' or genomes[ref].seq[(pos-1)] == '-':
                if (vcf['calldata/AD'][i,0,0])/(vcf['calldata/DP'][i,0]) < 0.5: # Checks that variant is a minority variant
                    mapping['position'].append(pos)
                    mapping['variant'].append(vcf['variants/ALT'][i,0])
                    mapping['coverage'].append(int(vcf['calldata/DP'][i,0]))
                    mapping['frequency'].append(float(vcf['calldata/AD'][i,0,0]/vcf['calldata/DP'][i,0]))
                elif (vcf['calldata/RD'][i,0])/(vcf['calldata/DP'][i,0]) >= 0.01: # Checks that reference variant meets cutoff
                    mapping['position'].append(pos)
                    mapping['variant'].append(vcf['variants/REF'][i])
                    mapping['coverage'].append(int(vcf['calldata/DP'][i,0]))
                    mapping['frequency'].append(float(vcf['calldata/RD'][i,0]/vcf['calldata/DP'][i,0]))

            else:
                if ((vcf['calldata/AD'][i,0,0])/(vcf['calldata/DP'][i,0]) >= 0.01 and 0.1 <= vcf['calldata/ADF'][i,0]/vcf['calldata/AD'][i,0,0] <= 0.9) or ((vcf['calldata/RD'][i,0])/(vcf['calldata/DP'][i,0]) >= 0.01 and 0.1 <= vcf['calldata/RDF'][i,0]/vcf['calldata/RD'][i,0] <= 0.9):
                    print('\nCheck ' + file + ' Genome does not match reference or variant.')
                    print('Strain: ' + ref)
                    print('Position: ' + str(pos))
                    print('Genome: '+ genomes[ref].seq[pos-1]).upper()
                    print('Ref: ' + vcf['variants/REF'][i])
                    print('Alt: ' + vcf['variants/ALT'][i,0])
    return mapping

def create_snvs(mfile, vcfs, genomes):
    '''
    Creates dictionary containing SNVs for each vcf file.
    '''
    with open(mfile, 'r') as f:
        metadata = pd.read_csv(f, sep = '\t')
    snvs = {}
    for index, row in metadata.iterrows():
        id = row['id']
        file = [path for path in vcfs if id in path]
        if len(file) > 0:
            reference = row['reference']
            snvs[id] = {}
            snvs[id] = check_variants(file[0], reference, genomes)
    return snvs

def add_total(snvs):
    '''
    Adds dictionary entry with total number of SNVs for each sample at each maf.
    '''
    for sample in snvs.keys():
        snvs[sample]['total'] = len(snvs[sample]['position'])
    return snvs

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Makes JSON of iSNVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True, help='tsv with sample names and ref')
    parser.add_argument('--sequences', type=str, required=True, help = 'fasta file of consensus genomes')
    parser.add_argument('--vcf', nargs = '+', type=str, required=True, help = 'directory containing vcf files')
    parser.add_argument('--output', type=str, required=True, help = 'location of output json')
    args = parser.parse_args()

    # Loads dictionary of consensus genomes
    genomes = load_genomes(args.sequences)
    #check_genomes(genomes, args.metadata)

    # Makes dictionary of iSNVs
    snvs = create_snvs(args.metadata, args.vcf, genomes)

    # Adds entry to iSNVs dictionary
    snvs_dict = add_total(snvs)

    # Writes iSNVs to JSON
    with open(args.output, 'wt') as f:
        json.dump(snvs_dict, f, indent=1)
