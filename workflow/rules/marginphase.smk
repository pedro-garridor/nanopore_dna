'''
Nanopore DNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

configfile: "config/config.yaml"

rule combine_vcfs:
    input:
        snvs=config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.phased.vcf.gz',
        svs=config['outdir']+'/sniffles2/{sample}.vcf.gz',# config['outdir']+'/HapDiff/{sample}'
    output:
        temp(config['outdir']+'/marginPhase/{sample}.vcf')
    conda:
        '../envs/samtools.yml'
    shell:
        '''
        bcftools concat \
            -a \
            {input.snvs} \
            {input.svs} \
            -o {output}
        '''

rule marginphase:
    input:
        bam=config['outdir']+'/BAM/{sample}.bam',
        bai=config['outdir']+'/BAM/{sample}.bam.bai',
        genome=config['genome_fa'],
        vcf=config['outdir']+'/marginPhase/{sample}.vcf',
        config='resources/allParams.phase_vcf.ont.sv.json'
    output:
        temp(directory(config['outdir']+'/marginPhase/{sample}'))
    threads:
        workflow.cores/4
    singularity:
        'docker://kishwars/pepper_deepvariant:r0.8'
    shell:
        '''
        margin phase \
            {input.bam} \
            {input.genome} \
            {input.vcf} \
            {input.config} \
            -t {threads} \
            -o {output} \
            -M
        '''

rule final_vcf:
    input:
        config['outdir']+'/marginPhase/{sample}'
    output:
        config['outdir']+'/VCF/{sample}.vcf.gz'
    conda:
        '../envs/samtools.yml'
    shell:
        'bgzip {input}'

rule final_vcf_idx:
    input:
        config['outdir']+'/VCF/{sample}.vcf.gz'
    output:
        config['outdir']+'/VCF/{sample}.vcf.gz.tbi'
    conda:
        '../envs/samtools.yml'
    shell:
        'tabix {input}'
