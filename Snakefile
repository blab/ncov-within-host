'''
To start, this pipeline calls minor allele variants from pileup files for SARS-CoV-2.
It will likely be expanded to do a lot more in future.
'''
import glob

SAMPLES = glob_wildcards('/fh/fast/bedford_t/seattleflu/assembly-ncov/{batch}/process/mpileup/sars-cov-2/{sample}.pileup').sample

rule all:
    input:
        snvs_df = 'results/snvs.json'

rule download:
    message: 'Downloading metadata and fasta files from S3'
    output:
        sequences = 'data/sequences_global.fasta',
        metadata = 'data/metadata_global.tsv'
    shell:
        '''
        aws s3 cp s3://nextstrain-ncov-private/metadata.tsv.gz - | gunzip -cq > {output.metadata}
        aws s3 cp s3://nextstrain-ncov-private/sequences.fasta.gz - | gunzip -cq > {output.sequences}
        '''

rule filter:
    message: 'Filtering to within-host Seattle Flu Study dataset'
    input:
        sequences = rules.download.output.sequences,
        metadata = rules.download.output.metadata
    output:
        sequences = 'results/genomes.fasta'
    params:
        min_length = 25000,
        exclude_where = 'submitting_lab!="Seattle Flu Study"',
        include_where = 'submitting_lab="Seattle Flu Study, University of Washington Medical Center"'
    shell:
        '''
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-length {params.min_length} \
            --exclude-where {params.exclude_where} \
            --include-where {params.include_where} \
            --output {output.sequences}
        '''

rule align:
    message: 'Aligning sequences to Wuhan/Hu-1/2019'
    input:
        sequences = rules.filter.output.sequences,
        reference = 'config/sars-cov-2.fasta'
    output:
        sequences = 'results/aligned.fasta'
    shell:
        '''
        augur align \
        --sequences {input.sequences} \
        --reference-sequence {input.reference} \
        --remove-reference \
        --output {output.sequences} \
        --remove-reference \
        --nthreads 4
        '''

rule clean_metabase:
    message: 'Cleaning metadata downloaded from Seattle Flu Metabase'
    input:
        metadata = 'data/metadata_metabase.tsv'
    output:
        metadata = 'data/metadata_metabase_clean.tsv'
    shell:
        '''
        awk 'BEGIN{{FS=OFS="\t"}}{{print $2,$4,$6,$7,$8,$11}}' {input.metadata} \
        | sed 's/avg_hcov19_crt/avg_ct/; s/"//g; s/\]//g; s/\[//g' \
        > {output.metadata}
        '''

def list_wadoh(wildcards):
    return glob.glob('data/wadoh_metadata/*.txt')

rule clean_wadoh:
    message: 'Cleaning and compiling WA-DoH metadata'
    input:
        metadata = list_wadoh
    output:
        metadata = 'data/metadata_wadoh_clean.tsv'
    shell:
        '''
        python scripts/clean_wadoh_metadata.py \
        --metadata {input.metadata} \
        --output {output.metadata}
        '''

rule concat_metadata:
    message: 'Combining Metabase, WA-DoH, and strains metadata'
    input:
        strains = '/fh/fast/bedford_t/seattleflu/aspera-data-backup/Flu/hcov19-batch-nwgc-id-strain.csv',
        wadoh = rules.clean_wadoh.output.metadata,
        global_metadata = rules.download.output.metadata,
        metabase = rules.clean_metabase.output.metadata
    output:
        metadata = 'results/metadata.tsv'
    shell:
        '''
        python scripts/concat_metadata.py \
        --strains {input.strains} \
        --wadoh {input.wadoh} \
        --global-metadata {input.global_metadata} \
        --metabase {input.metabase} \
        --output {output.metadata}
        '''

def list_pileups(wildcards):
    return glob.glob('/fh/fast/bedford_t/seattleflu/assembly-ncov/*/process/mpileup/sars-cov-2/' + wildcards.sample +'.pileup')

rule call_snvs:
    message: 'Calling SNVs compared to Wuhan-1 using varscan'
    input:
        pileup = list_pileups
    output:
        vcf = 'results/vcf/{sample}.vcf'
    params:
        min_cov = 100,
        phred = 30,
        min_var_freq = 0.01
    shell:
        '''
        varscan mpileup2snp \
        {input.pileup} \
        --min-coverage {params.min_cov} \
        --min-avg-qual {params.phred} \
        --min-var-freq {params.min_var_freq} \
        --strand-filter 1 \
        --output-vcf 1 \
        > {output.vcf}
        '''

rule validate_snvs:
    message: 'Generating iSNVs for each sample'
    input:
        metadata = rules.concat_metadata.output.metadata,
        sequences = rules.align.output.sequences,
        vcfs = expand('results/vcf/{sample}.vcf', sample=SAMPLES)
    output:
        snvs = 'results/snvs.json'
    shell:
        '''
        python scripts/validate_snvs.py \
        --metadata {input.metadata} \
        --sequences {input.sequences} \
        --vcf {input.vcfs} \
        --output {output.snvs}
        '''

#rule plot_snvs_ct:
#    message: 'Plotting iSNVs vs. Ct'
#    input:
#        snvs = rules.validate_snvs.output.snvs,
#        metadata = rules.concat_metadata.output.metadata
#    output:
#        plot = 'figures/ct-snvs.pdf'
#    shell:
#        '''
#        python scripts/plot_snvs_ct.py \
#        --snvs {input.snvs} \
#        --metadata {input.metadata} \
#        --output {output.plot}
#        '''

#rule plot_freq_ct:
#    message: 'Plotting frequency of iSNVs for SCAN samples colored by Ct'
#    input:
#        snvs = rules.validate_snvs.output.snvs,
#        metadata = rules.concat_metadata.output.metadata
#    output:
#        plot = 'figures/ct-freq-SCAN.pdf'
#    shell:
#        '''
#        python scripts/plot_freq_ct.py \
#        --snvs {input.snvs} \
#        --metadata {input.metadata} \
#        --output {output.plot}
#        '''

#rule plot_snvs_freq_ct:
#    message: 'Plotting total and frequency of iSNVs vs. Ct for SCAN samples'
#    input:
#        snvs = rules.validate_snvs.output.snvs,
#        metadata = rules.concat_metadata.output.metadata
#    output:
#        plot = 'figures/ct-snvs-freq-SCAN.pdf'
#    shell:
#        '''
#        python scripts/plot_snvs_freq_ct.py \
#        --snvs {input.snvs} \
#        --metadata {input.metadata} \
#        --output {output.plot}
#        '''

rule choose_duplicates:
    message: 'Choosing samples to sequence in duplicate'
    input:
        metadata = rules.concat_metadata.output.metadata,
        snvs = rules.validate_snvs.output.snvs
    output:
        duplicates = 'results/to_duplicate.tsv'
    shell:
        '''
        python scripts/choose_duplicates.py \
        --metadata {input.metadata} \
        --snvs {input.snvs} \
        --output {output.duplicates}
        '''

#rule create_snvs_df:
#    message: 'Saves DF containing all SNVs specified by params.origin & params.ct_cutoff'
#    input:
#        metadata = 'data/metadata.tsv',
#        snvs = 'results/snvs.json'
#    output:
#        snvs = 'results/snvs_sfs.tsv'
#    params:
#        origin = 'sfs',
#        ct_cutoff = 22
#    shell:
#        '''
#        python scripts/create_snvs_df.py \
#        --metadata {input.metadata} \
#        --snvs {input.snvs} \
#        --origin {params.origin} \
#        --ct-cutoff {params.ct_cutoff} \
#        --output {output.snvs}
#        '''
