'''
Nanopore DNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

rule shasta:
    input:
        config['outdir']+'/FASTQ/{sample}.fq'
    output:
        temp(directory(config['outdir']+'/SHASTA/{sample}'))
    params:
        config='resources/Nanopore-R10-Fast-Nov2022.conf'
    threads:
        workflow.cores/4
    conda:
        './envs/shasta.yml'
    shell:
        '''
        shasta assemble\
            --input {input} \
            --assemblyDirectory {output} \
            --config {params.config} \
            --threads {threads}
        '''

rule minimap2_diploid:
    input:
        diploid_ref=config['outdir']+'/SHASTA/{sample}',
        fq=config['outdir']+'/FASTQ/{sample}.fq'
    output:
        temp(config['outdir']+'/HapDup/{sample}/{sample}.sam')
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
            {input.diploid_ref}/Assembly.fasta\
            {input.fq} > {output}
        '''

rule samtools_sort_diploid:
    input:
        config['outdir']+'/HapDup/{sample}/{sample}.sam'
    output:
        temp(config['outdir']+'/HapDup/{sample}/{sample}.bam')
    conda:
        '../envs/samtools.yml'
    shell:  
        'samtools sort {input} -@ {threads} -o {output}'

rule samtools_index_diploid:
    input: 
        config['outdir']+'/HapDup/{sample}/{sample}.bam'
    output:
        temp(config['outdir']+'/HapDup/{sample}/{sample}.bam.bai')
    conda:
        '../envs/samtools.yml'
    shell:  
        'samtools index {input}'

rule hapdup:
    input:
        diploid_ref=config['outdir']+'/SHASTA/{sample}',
        bam=config['outdir']+'/HapDup/{sample}/{sample}.bam'
        bai=config['outdir']+'/HapDup/{sample}/{sample}.bam.bai'
    output:
        temp(directory(config['outdir']+'/HapDup/{sample}'))
    threads:
        workflow.cores/4
    singularity:
        'docker://mkolmogo/hapdup'
    shell:
        '''
        hapdup \
            --assembly {input.diploid_ref}/Assembly.fasta \
            --bam {input.bam} \
            --out-dir {output} \
            -t {threads}
        '''
