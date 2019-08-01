#!/usr/bin/env python3

#map ogbf numbers to samples names so that in future can refer to sample names

def resolve_input_reads(wildcards):
  sample_number = sample_name_to_ogbf_number[wildcards.sample]
  return({
  'L1R1': ('data/reads/H2NVLBCX3-4261-'
           '{sample_number:02}-01-01_'
           'S{sample_number}_L001_R1_001.fastq.gz'),
  'L1R2': ('data/reads/H2NVLBCX3-4261-'
           '{sample_number:02}-01-01_'
           'S{sample_number}_L001_R2_001.fastq.gz'),
  'L2R1': ('data/reads/H2NVLBCX3-4261-'
           '{sample_number:02}-01-01_'
           'S{sample_number}_L002_R1_001.fastq.gz'),
  'L2R2': ('data/reads/H2NVLBCX3-4261-'
           '{sample_number:02}-01-01_'
           'S{sample_number}_L002_R2_001.fastq.gz')})

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

rule target:
  input:
    expand('output/010_trimmed/{sample}_R1.fq.gz', 
           sample=list(sample_name_to_ogbf_number.keys()))

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
