#!/usr/bin/env python3

import multiprocessing

#map ogbf numbers to samples names so that in future can refer to sample names

def resolve_input_reads(wildcards):
  barcode_number = sample_name_to_ogbf_number[wildcards.sample]
  library_number = barcode_number if barcode_number <= 5 else barcode_number - 1
  return({
  'L1R1': ('data/reads/H2NVLBCX3-4261-'
           f'{barcode_number:02}-01-01_'
           f'S{library_number}_L001_R1_001.fastq.gz'),
  'L1R2': ('data/reads/H2NVLBCX3-4261-'
           f'{barcode_number:02}-01-01_'
           f'S{library_number}_L001_R2_001.fastq.gz'),
  'L2R1': ('data/reads/H2NVLBCX3-4261-'
           f'{barcode_number:02}-01-01_'
           f'S{library_number}_L002_R1_001.fastq.gz'),
  'L2R2': ('data/reads/H2NVLBCX3-4261-'
           f'{barcode_number:02}-01-01_'
           f'S{library_number}_L002_R2_001.fastq.gz')})

sample_name_to_ogbf_number = {
  'AI': 1,
  'AI_b': 2,
  'H1_BF': 3,
  'H1_BYc': 4,
  'H3_BF': 5,
  'H2_OxFc': 7,
  'H2_OxY': 8,
  'H3_A2': 9,
  'H3_N3': 10,
  'H5_A4': 11,
  'H5_N3': 12,
  'H1': 13,
  'H1b': 14,
  'H2b': 15,
  'H3b': 16}

bbmap = 'shub://TomHarrop/singularity-containers:bbmap_38.50b'
salmon = 'local_containers/salmon_0.14.1.sif'
salmontools = 'local_containers/salmontools_23eac84.sif'
bioconductor = 'local_containers/bioconductor_3.9.sif'

rule target:
  input:
    expand('output/030_salmon/{sample}/quant.sf', 
           sample=list(sample_name_to_ogbf_number.keys()))

rule Generate_DESeq_object:
  input:
    quant_files=expand('output/030_salmon/{sample}/quant.sf',
                      sample=list(sample_name_to_ogbf_number.keys())),
     gff = 'data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.gff',
     sample_summary = 'data/sample_summary.csv'
  output:
    dds = 'output/040_DESeq/dds.Rds'
  log:
    'output/logs/Generate_DESeq_object.log'
  script:
    'source/Generate_DESeq_object.R'
  singularity:
    bioconductor

rule quantification:
  input:
    R1='output/010_trimmed/{sample}_R1.fq.gz',
    R2='output/010_trimmed/{sample}_R2.fq.gz',
    index='output/005_index'
  output:
    'output/030_salmon/{sample}/quant.sf'
  params:
    outdir='output/030_salmon/{sample}'
  log:
    'output/logs/030_salmon/{sample}.log'
  threads:
    multiprocessing.cpu_count()
  singularity:
    salmon
  shell:
    'salmon quant '
    '--libType ISR '
    '--index {input.index} '
    '--mates1 {input.R1} '
    '--mates2 {input.R2} '
    '--output {params.outdir} '
    '--threads {threads} '
    '--validateMappings '
    '--gcBias '
    '&> {log} '

rule trim_adaptors:
  input:
    R1='output/000_merged/{sample}_R1.fq.gz',
    R2='output/000_merged/{sample}_R2.fq.gz'
  output:
    R1='output/010_trimmed/{sample}_R1.fq.gz',
    R2='output/010_trimmed/{sample}_R2.fq.gz'
  params:
    adapters='/adapters.fa'
  log:
    r='output/logs/010_trim/{sample}_repair.log',
    b='output/logs/010_trim/{sample}_bbduk.log'
  threads:
    1
  singularity:
    bbmap
  shell:
    'repair.sh '
    'in={input.R1} '
    'in2={input.R2} '
    'out=stdout.fastq '
    '2>{log.r} '
    '| '
    'bbduk.sh '
    'threads={threads} '
    'in=stdin.fastq '
    'int=t '
    'out={output.R1} '
    'out2={output.R2} '
    'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
    'ref={params.adapters} '
    '2>{log.b} '

rule merge_lanes:
  input:
    unpack(resolve_input_reads)
  output:
    R1=temp('output/000_merged/{sample}_R1.fq.gz'),
    R2=temp('output/000_merged/{sample}_R2.fq.gz')
  singularity:
    bbmap
  shell:
    'cat {input.L1R1} {input.L2R1} > {output.R1} & '
    'cat {input.L1R2} {input.L2R2} > {output.R2} & '
    'wait'

rule generate_index:
    input:
        transcriptome = 'output/020_ref/gentrome.fa',
        decoys = 'output/020_ref/decoys.txt'
    output:
        directory('output/005_index')
    log:
        'output/logs/005_index/generate_index.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        salmon
    shell:
        'salmon index '
        '--transcripts {input.transcriptome} '
        '--index {output} '
        '--threads {threads} '
        '--decoys {input.decoys} '
        '&> {log}'

rule generate_decoy_trancriptome:
    input:
        fasta = 'data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.fna',
        transcriptome = 'data/ref/GCF_003254395.2_Amel_HAv3.1_rna.fna',
        gff = 'data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.gff'
    output:
        'output/020_ref/gentrome.fa',
        'output/020_ref/decoys.txt'
    params:
        outdir = 'output/020_ref'
    log:
        'output/logs/020_ref/generate_decoy_trancriptome.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        salmontools
    shell:
        'generateDecoyTranscriptome.sh '
        '-j {threads} '
        '-b /usr/bin/bedtools '
        '-m /usr/local/bin/mashmap '
        '-a {input.gff} '
        '-g {input.fasta} '
        '-t {input.transcriptome} '
        '-o {params.outdir} '
        '&> {log}'
