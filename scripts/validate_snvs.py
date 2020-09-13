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

def strains_to_samples(tsv, fasta):
    '''
    Renames keys of dictionary of consensus genomes with nwgc_id.
    Drops genomes without nwgc_id.
    '''
    with open(tsv) as tfile:
        metadata = pd.read_csv(tfile, sep = '\t')

    genomes = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    strains = [strain for strain in genomes.keys()]
    keep = list(set(strains) & set(metadata['strain']))
    drop = list(set(strains).difference(metadata['strain']))

    for strain in keep:
        nwgc_id = metadata.loc[metadata.strain == strain, 'nwgc_id'].values
        if len(nwgc_id) > 1: # Assumes that there is only duplicate sequencing of samples. If more, e.g. triplicate, modify code.
            genomes[nwgc_id[0]] = genomes[strain]
            genomes[nwgc_id[1]] = genomes.pop(strain)
        else:
            genomes[nwgc_id[0]] = genomes.pop(strain)

    for strain in drop:
        genomes.pop(strain)
        print(strain + " NOT in " + tsv)
    return genomes

def check_variants(file, nwgc_id, genomes):
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

            if genomes[nwgc_id].seq[(pos-1)] == vcf['variants/REF'][i]: # Minus 1 corrects for position numbering beginning at 0 for SeqIO records.
                if 0.1 <= (vcf['calldata/ADF'][i,0])/(vcf['calldata/AD'][i,0,0]) <= 0.9: # Checks that variant is supported by both strands
                    mapping['position'].append(pos)
                    mapping['variant'].append(vcf['variants/ALT'][i,0])
                    mapping['coverage'].append(int(vcf['calldata/DP'][i,0]))
                    mapping['frequency'].append(float(vcf['calldata/AD'][i,0,0]/vcf['calldata/DP'][i,0]))

            elif genomes[nwgc_id].seq[(pos-1)] == vcf['variants/ALT'][i,0]: # Minus 1 corrects for position numbering beginning at 0 for SeqIO records.
                if (vcf['calldata/RD'][i,0])/(vcf['calldata/DP'][i,0]) >= 0.01: # Checks that reference variant meets cutoff
                    if 0.1 <= (vcf['calldata/RDF'][i,0])/(vcf['calldata/RD'][i,0]) <= 0.9: # Checks that variant is supported by both strands
                        mapping['position'].append(pos)
                        mapping['variant'].append(vcf['variants/REF'][i])
                        mapping['coverage'].append(int(vcf['calldata/DP'][i,0]))
                        mapping['frequency'].append(float(vcf['calldata/RD'][i,0]/vcf['calldata/DP'][i,0]))

            elif genomes[nwgc_id].seq[(pos-1)] == 'N':
                if (vcf['calldata/AD'][i,0,0])/(vcf['calldata/DP'][i,0]) < 0.5: # Checks that variant is a minority variant
                    if 0.1 <= (vcf['calldata/ADF'][i,0])/(vcf['calldata/AD'][i,0,0]) <= 0.9: # Checks that variant is supported by both strands
                        mapping['position'].append(pos)
                        mapping['variant'].append(vcf['variants/ALT'][i,0])
                        mapping['coverage'].append(int(vcf['calldata/DP'][i,0]))
                        mapping['frequency'].append(float(vcf['calldata/AD'][i,0,0]/vcf['calldata/DP'][i,0]))
                elif (vcf['calldata/RD'][i,0])/(vcf['calldata/DP'][i,0]) >= 0.01: # Checks that reference variant meets cutoff
                    if 0.1 <= (vcf['calldata/RDF'][i,0])/(vcf['calldata/RD'][i,0]) <= 0.9: # Checks that variant is supported by both strands
                            mapping['position'].append(pos)
                            mapping['variant'].append(vcf['variants/REF'][i])
                            mapping['coverage'].append(int(vcf['calldata/DP'][i,0]))
                            mapping['frequency'].append(float(vcf['calldata/RD'][i,0]/vcf['calldata/DP'][i,0]))

            else:
                if ((vcf['calldata/AD'][i,0,0])/(vcf['calldata/DP'][i,0]) >= 0.01 and 0.1 <= vcf['calldata/ADF'][i,0]/vcf['calldata/AD'][i,0,0] <= 0.9) or ((vcf['calldata/RD'][i,0])/(vcf['calldata/DP'][i,0]) >= 0.01 and 0.1 <= vcf['calldata/RDF'][i,0]/vcf['calldata/RD'][i,0] <= 0.9):
                    print('\nCheck ' + str(nwgc_id) + '.vcf. Genome does not match reference or variant.')
                    print('Position: ' + str(pos))
                    print('Genome: '+ genomes[nwgc_id].seq[pos-1])
                    print('Ref: ' + vcf['variants/REF'][i])
                    print('Alt: ' + vcf['variants/ALT'][i,0])
    return mapping

def create_snvs(vcfs, genomes):
    '''
    Creates dictionary containing SNVs for each vcf file.
    '''
    snvs = {}
    for file in vcfs:
        path_list = file.split('/')
        nwgc_id = int(path_list[2].split('.')[0])

        if nwgc_id in genomes.keys():
            snvs[nwgc_id] = {}
            snvs[nwgc_id] = check_variants(file, nwgc_id, genomes)
    #    else:
    #        print('For ' + file + ', sample is not in fasta.')
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

    parser.add_argument('--metadata', type=str, required=True, help='tsv with sample names and nwgc_id')
    parser.add_argument('--sequences', type=str, required=True, help = 'fasta file of consensus genomes')
    parser.add_argument('--vcf', nargs = '+', type=str, required=True, help = 'directory containing vcf files')
    parser.add_argument('--output', type=str, required=True, help = 'location of output json')
    args = parser.parse_args()

    # Loads dictionary of consensus genomes
    genomes = strains_to_samples(args.metadata, args.sequences)

    # Makes dictionary of iSNVs
    snvs = create_snvs(args.vcf, genomes)

    # Adds entry to iSNVs dictionary
    snvs_dict = add_total(snvs)

    # Writes iSNVs to JSON
    with open(args.output, 'wt') as f:
        json.dump(snvs_dict, f, indent=1)
