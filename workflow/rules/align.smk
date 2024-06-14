'''
Nanopore DNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

rule fastq:
    input:
        config['input_dir']+'/{sample}.bam'
    output:
        temp(config['outdir']+'/FASTQ/{sample}.fq')
    conda: '../envs/samtools.yml'
    shell:
        '''
        samtools fastq \
            -T* \
            {input} > {output}
        '''
# -TMm,Ml,MM,ML

rule minimap2:
    input:
        genome=config['genome_dir']+'/hg38.fa',
        fq=config['outdir']+'/FASTQ/{sample}.fq'
    output:
        temp(config['outdir']+'/BAM/{sample}.sam')
    threads:
        workflow.cores/2
    conda:
        '../envs/minimap2.yml'
    shell:
        '''
        minimap2 \
            -ax map-ont \
            --MD \
            -k 14 \
            -t {threads} \
            -y \
            {input} > {output}
        '''

rule samtools_sort:
    input:
        config['outdir']+'/BAM/{sample}.sam'
    output:
        protected(config['outdir']+'/BAM/{sample}.bam')
    threads:
        workflow.cores/4
    conda:
        '../envs/samtools.yml'
    shell:
        'samtools sort {input} -@ {threads} -o {output}'

rule samtools_index:
    input:
        config['outdir']+'/BAM/{sample}.bam',
    output:
        protected(config['outdir']+'/BAM/{sample}.bam.bai')
    conda:
        '../envs/samtools.yml'
    shell:
        'samtools index {input}'
