'''
Constructs df with all SNVs annotated with amino acid changes from samples in metadata tsv.
'''
import argparse
import json
import pandas as pd
from Bio import SeqIO
import math
from Bio.Seq import Seq

def load_metadata(file):
    '''
    Loads metadata tsv as df.
    '''
    with open(file) as tfile:
        df = pd.read_csv(tfile, sep = '\t')
        df['nwgc_id'] = df.nwgc_id.astype('str')
    return df[['strain', 'id', 'nwgc_id', 'date', 'avg_ct', 'primers', 'symptom_onset', 'address_identifier', 'individual_identifier', 'sample_identifier', 'age_bin', 'age', 'location', 'puma', 'origin', 'batch']]

def load_snvs(file):
    '''
    Loads SNVs dictionary.
    '''
    with open(file) as jfile:
        snvs = json.load(jfile)
    return snvs

def load_genomes(file):
    '''
    Loads fasta as dictionary.
    '''
    with open(file) as fasta:
        genomes = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    return genomes

def load_gb(file):
    '''
    Loads Genbank file.
    '''
    with open(file) as handle:
        record = SeqIO.read(handle, "genbank")
    return record

def create_snvs_df(file, metadata, snvs):
    '''
    Constructs df with all SNVs
    '''
    with open(file) as tfile:
        samples_df = pd.read_csv(tfile, sep = '\t')
        samples_df['nwgc_id'] = samples_df.nwgc_id.astype('str')
    samples = samples_df.nwgc_id.to_list()

    positions = []
    alts = []
    frequencies = []
    coverage = []
    ids = []
    for sample in samples:
        positions.extend(snvs[sample]['position'])
        alts.extend(snvs[sample]['variant'])
        frequencies.extend(snvs[sample]['frequency'])
        coverage.extend(snvs[sample]['coverage'])
        ids.extend([sample]*len(snvs[sample]['position']))
    df = pd.DataFrame()
    df['nwgc_id'] = ids
    df['position'] = positions
    df['alt'] = alts
    df['coverage'] = coverage
    df['frequency'] = frequencies
    df = df.merge(metadata, on='nwgc_id', how='left')
    return df

def add_consensus_bp(df, genomes):
    '''
    Adds nucleotide in consensus genome at variant location to df.
    '''
    consensus = []
    for sample, pos in zip(df['strain'], df['position']):
        bp = genomes[sample].seq[(pos-1)].upper()
        consensus.append(bp)
    df['consensus'] = consensus
    return df

def annotate_variants(df, record, genomes):
    '''
    Adds amino acid changes to variant.
    '''
    protein = []
    aa_position = []
    aa_consensus = []
    aa_alt = []

    for i in df.index:
        pos = int(df.loc[i, 'position'])
        strain = df.loc[i, 'strain']
        alt = df.loc[i, 'alt']
        proteins = []
        aa_positions = []
        aa_cns = []
        aa_alts = []

        for feature in record.features:
            f = feature.location
            if feature.type != 'gene' and feature.type != 'source' and pos in feature:
                if feature.type == 'CDS':
                    proteins.append(feature.qualifiers['gene'][0])
                    aa_pos = math.ceil((pos - f.start)/3)
                    r = (pos - f.start) % 3
                    if r == 0:
                        codon = genomes[strain].seq[(pos-3):pos]
                        variant_codon = Seq(str(codon[0:2] + alt))
                    else:
                        codon = genomes[strain].seq[(pos-r):(pos+(3-r))]
                        if r == 1:
                            variant_codon = Seq(alt + str(codon[1:3]))
                        elif r == 2:
                            variant_codon = Seq(str(codon[0]) + alt + str(codon[2]))
                    aa = str(codon.translate())
                    aa_var = str(variant_codon.translate())
                    aa_positions.append(aa_pos)
                    aa_cns.append(aa)
                    aa_alts.append(aa_var)
                else:
                    proteins.append(feature.type)

        if len(proteins) == 0:
            proteins.append('non-coding')
        str_proteins = ";".join(proteins)
        str_aa_pos = [str(pos) for pos in aa_positions]
        str_positions = ";".join(str_aa_pos)
        str_aa = ";".join(aa_cns)
        str_var = ";".join(aa_alts)
        protein.append(str_proteins)
        aa_position.append(str_positions)
        aa_consensus.append(str_aa)
        aa_alt.append(str_var)

    df['protein'] = protein
    df['aa_position'] = aa_position
    df['aa_consensus'] = aa_consensus
    df['aa_alt'] = aa_alt
    return df

def add_mutation_type(df):
    '''
    Adds mutation type to df
    '''
    mutations = []
    for row in df.index:
        protein = df.loc[row, 'protein']
        if protein == "5'UTR" or protein == "3'UTR" or protein == 'non-coding':
            mutations.append('NC')
        elif df.loc[row, 'aa_consensus'] == df.loc[row, 'aa_alt']:
            mutations.append('S')
        else:
            mutations.append('NS')
    df['mut_type'] = mutations
    return df

def add_variant(df):
    '''
    Adds nt variant column derived from consensus, position, and alt.
    '''
    df['variant'] = df['consensus'].str.cat(df['position'].astype('str'),sep="").str.cat(df['alt'],sep="")
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

    parser.add_argument('--metadata', type=str, required=True, help='tsv file containing metadata')
    parser.add_argument('--snvs', type=str, required=True, help='iSNVs json')
    parser.add_argument('--sequences', type=str, required=True, help='Fasta file with aligend genomes')
    parser.add_argument('--reference', type=str, required=True, help='SARS-CoV-2 reference genbank file')
    parser.add_argument('--ids', type=str, required=True, help='tsv with list of nwgc ids')
    parser.add_argument('--output', type=str, required=True, help='location of output tsv')
    args = parser.parse_args()

    # Load files
    metadata = load_metadata(args.metadata)
    snvs = load_snvs(args.snvs)
    genomes =  load_genomes(args.sequences)
    genbank = load_gb(args.reference)

    # Create SNVS df
    snvs_df = create_snvs_df(args.ids, metadata, snvs)

    # Add to df
    snvs_df = add_consensus_bp(snvs_df, genomes)
    snvs_df = annotate_variants(snvs_df, genbank, genomes)
    snvs_df = add_mutation_type(snvs_df)
    snvs_df = add_variant(snvs_df)

    # Reorder columns
    columns_order = ['strain', 'nwgc_id', 'variant', 'protein', 'mut_type', 'aa_consensus', 'aa_position', 'aa_alt', 'frequency', 'coverage', 'avg_ct', 'date', 'symptom_onset', 'address_identifier']
    snvs_df = reorder_cols(snvs_df, columns_order)

    # Write out tsv
    with open(args.output, 'w') as f:
        snvs_df.to_csv(f, sep = '\t', index=False)
