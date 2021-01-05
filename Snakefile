'''
To start, this pipeline calls minor allele variants from pileup files for SARS-CoV-2.
It will likely be expanded to do a lot more in future.
'''
configfile: "config/config.yaml"

import glob

rule all:
    input:
        snvs_df = 'results/snvs.json'

rule download_metadata:
    message: "Downloading GISAID SARS-CoV-2 metadata from S3"
    output:
        metadata = "data/metadata_gisaid.tsv"
    shell:
        """
        aws s3 cp s3://nextstrain-ncov-private/metadata.tsv.gz - | gunzip -cq >{output.metadata:q}
        """

rule download_aligned:
    message: "Downloading aligned SARS-CoV-2 GISAID fasta files from S3 bucket"
    output:
        sequences = "data/aligned_gisaid.fasta"
    params:
        compression = "xz",
        deflate = "xz -dcq"
    shell:
        """
        aws s3 cp s3://nextstrain-ncov-private/aligned.fasta.{params.compression} - | {params.deflate} > {output.sequences:q}
        """

rule filter:
    message: 'Filtering to within-host Seattle Flu Study dataset'
    input:
        sequences = rules.download_aligned.output.sequences,
        metadata = rules.download_metadata.output.metadata
    output:
        sequences = 'results/genomes.fasta'
    params:
        exclude_where = 'submitting_lab!="Seattle Flu Study"',
        include_where = 'submitting_lab="Seattle Flu Study, University of Washington Medical Center"'
    shell:
        '''
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude-where {params.exclude_where} \
            --include-where {params.include_where} \
            --output {output.sequences}
        '''

rule clean_metabase:
    message: 'Cleaning metadata downloaded from Seattle Flu Metabase'
    input:
        metadata = 'data/metadata_metabase.tsv'
    output:
        metadata = 'data/metadata_metabase_clean.tsv'
    shell:
        '''
        awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$3,$5,$6,$8,$9,$10,$13}}' {input.metadata} \
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
    message: 'Combining Metabase, WA-DoH, and strains metadata & filtering to genomes'
    input:
        strains = '/fh/fast/bedford_t/seattleflu/aspera-data-backup/Flu/hcov19-batch-nwgc-id-strain.csv',
        wadoh = rules.clean_wadoh.output.metadata,
        global_metadata = rules.download_metadata.output.metadata,
        metabase = rules.clean_metabase.output.metadata,
        genomes = rules.filter.output.sequences
    output:
        metadata = 'results/metadata.tsv'
    shell:
        '''
        python scripts/concat_metadata.py \
        --strains {input.strains} \
        --wadoh {input.wadoh} \
        --global-metadata {input.global_metadata} \
        --metabase {input.metabase} \
        --genomes {input.genomes} \
        --output {output.metadata}
        '''

def list_pileups(wildcards):
    return glob.glob('/fh/scratch/delete10/bedford_t/ncov_pileups/' + wildcards.sample +'.pileup')

rule call_snvs:
    message: 'Calling SNVs compared to Wuhan-1 using varscan'
    input:
        pileup = list_pileups
    output:
        vcf = 'results/vcf_snvs/{sample}.vcf'
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
        --strand-filter 0 \
        --variants 1 \
        --output-vcf 1 \
        > {output.vcf}
        '''

rule call_indels:
    message: 'Calling SNVs compared to Wuhan-1 using varscan'
    input:
        pileup = list_pileups
    output:
        vcf = 'results/vcf_indels/{sample}.vcf'
    params:
        min_cov = 100,
        phred = 30,
        min_var_freq = 0.01
    shell:
        '''
        varscan mpileup2indel \
        {input.pileup} \
        --min-coverage {params.min_cov} \
        --min-avg-qual {params.phred} \
        --min-var-freq {params.min_var_freq} \
        --strand-filter 0 \
        --variants 1 \
        --output-vcf 1 \
        > {output.vcf}
        '''

rule get_indels:
    message: 'dummy rule to force indel calling'
    input:
        vcfs = expand('results/vcf_indels/{sample}.vcf', sample=config['samples'])
    output:
        text = 'results/indels_done.txt'
    shell:
        '''
        echo "Yeehaw" > {output.txt}
        '''

rule validate_snvs:
    message: 'Generating iSNVs for each sample'
    input:
        metadata = rules.concat_metadata.output.metadata,
        sequences = rules.filter.output.sequences,
        vcfs = expand('results/vcf_snvs/{sample}.vcf', sample=config['samples'])
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

rule construct_snvs_df:
    message: 'Creates df with all SNVs + annotations for samples in given tsv'
    input:
        metadata = 'results/metadata_{origin}_{dataset}.tsv',
        snvs = rules.validate_snvs.output.snvs,
        sequences = rules.filter.output.sequences,
        reference = 'config/reference_seq.gb'
    output:
        snvs = 'results/snvs_{origin}_{dataset}.tsv'
    shell:
        '''
        python scripts/create_snvs_df.py \
        --metadata {input.metadata} \
        --snvs {input.snvs} \
        --sequences {input.sequences} \
        --reference {input.reference} \
        --output {output.snvs}
        '''
