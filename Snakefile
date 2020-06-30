'''
To start, this pipeline calls minor allele variants from pileup files for SARS-CoV-2.
It will likely be expanded to do a lot more in future.
'''
import glob

configfile: 'config/config.yaml'

rule all:
    input:
        isnvs = 'results/snvs.json'

rule download:
    message: 'Downloading metadata and fasta files from S3'
    output:
        sequences = 'data/preliminary/sequences_global.fasta',
        metadata = 'data/preliminary/metadata_global.tsv'
    shell:
        '''
        aws s3 cp s3://nextstrain-ncov-private/metadata.tsv.gz - | gunzip -cq >{output.metadata}
        aws s3 cp s3://nextstrain-ncov-private/sequences.fasta.gz - | gunzip -cq > {output.sequences}
        '''

rule filter:
    message: 'Filtering to within-host Seattle Flu Study dataset'
    input:
        sequences = rules.download.output.sequences,
        metadata = rules.download.output.metadata
    output:
        sequences = 'data/genomes.fasta'
    params:
        min_length = 25000,
        exclude_where = 'submitting_lab!="Seattle Flu Study"'
    shell:
        '''
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-length {params.min_length} \
            --exclude-where {params.exclude_where} \
            --output {output.sequences}
        '''

#rule clean_metabase: Will add in when I sort through all the rest of concating metadata rules
#    message: 'Cleaning metadata downloaded from Seattle Flu Metabase'
#    input:
#        metadata = 'data/preliminary/metadata_metabase.tsv'
#    output:
#        metadata = 'data/preliminary/metadata_metabase_clean.tsv'
#    shell:
#        '''
#        awk 'BEGIN{{FS=OFS="\t"}}{{split($3,a,"T"); print $2, a[1],$4,$6,$7,$8,$11}}' {input.metadata} \
#        | sed 's/avg_hcov19_crt/ct/; s/"//g' \
#        > {output.metadata}
#        '''

def list_pileups(wildcards):
    return glob.glob('/fh/fast/bedford_t/seattleflu/assembly-ncov/*/process/mpileup/sars-cov-2/' + wildcards.sample +'.pileup')

rule call_snvs:
    message: 'Calling SNVs compared to Wuhan-1 using varscan'
    input:
        pileup = list_pileups
    output:
        vcf = 'data/vcf/maf-{freq}/{sample}.vcf'
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
        strains = 'data/sars-cov-2_batch_nwgc-id_strain.tsv',
        sequences = rules.filter.output.sequences,
        vcfs = expand('data/vcf/{sample}.vcf', sample=config['samples'])
    output:
        snvs = 'results/snvs.json'
    shell:
        '''
        python scripts/validate_snvs.py \
        --strain-id {input.strains} \
        --sequences {input.sequences} \
        --vcf {input.vcfs} \
        --output {output.snvs}
        '''

rule plot_snvs_ct:
    message: 'Plotting iSNVs vs. Ct'
    input:
        snvs = rules.validate_snvs.output.snvs,
        metadata = 'data/SCAN_within-host-metadata_2020-06-13.tsv'
    output:
        plot = 'figures/ct-snvs.pdf'
    shell:
        '''
        python scripts/plot_snvs_ct.py \
        --snvs {input.snvs} \
        --metadata {input.metadata} \
        --output {output.plot}
        '''
