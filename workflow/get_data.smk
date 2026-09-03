import os

configfile: 'config/config.yaml'

REF_DIR = os.path.join(config['path']['data_dir'], config['path']['reference_dir'], f'{SPECIES}_{BUILD}_{RELEASE}')
SPECIES = config['run']['ref_species']
BUILD = config['run']['ref_build']
RELEASE = config['run']['ref_release']

rule get_cdna:
    output: 
        f'{REF_DIR}/{SPECIES}.{BUILD}.cdna.all.fa.gz'
    shell:
        '''
        wget -P {REF_DIR} 'http://ftp.ensembl.org/pub/release-{RELEASE}/fasta/{SPECIES}/cdna/{SPECIES}.{BUILD}.cdna.all.fa.gz'
        '''

rule get_genome:
    output: 
        f'{REF_DIR}/{SPECIES}.{BUILD}.dna.toplevel.fa.gz'
        f'{REF_DIR}/{SPECIES}.{BUILD}.dna.toplevel.fa.gz.fai'
    shell:
        '''
        wget -P {REF_DIR} 'http://ftp.ensembl.org/pub/release-{RELEASE}/fasta/{SPECIES}/dna_index/{SPECIES}.{BUILD}.dna.toplevel.fa.gz'
        wget -P {REF_DIR}'http://ftp.ensembl.org/pub/release-{RELEASE}/fasta/{SPECIES}/dna_index/{SPECIES}.{BUILD}.dna.toplevel.fa.gz.fai'
        '''

rule get_annotation:
    output: 
        f'{REF_DIR}/{SPECIES}.{BUILD}.{RELEASE}.gtf.gz'
    shell:
        '''
        wget -P {REF_DIR} 'http://ftp.ensembl.org/pub/release-{RELEASE}/gtf/{SPECIES}/{SPECIES}.{BUILD}.{RELEASE}.gtf.gz'
        '''