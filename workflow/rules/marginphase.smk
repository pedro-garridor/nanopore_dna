'''
Nanopore DNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

rule combine_vcfs:
    input:
        snvs=config['outdir']+'/PMDV/{sample}',
        svs=config['outdir']+'/HapDiff/{sample}'
    output:
        temp(config['outdir']+'/marginPhase/{sample}.vcf')
    conda:
        '../envs/bcftools.yml'
    shell:
        '''
        bcftools concat \
            -a \
            {input.snvs}/{sample}_PMDV_FINAL.phased.vcf.gz \
            {input.svs}/hapdiff_unphased.vcf.gz -o {output}
        '''

rule marginphase:
    input:
        bam=config['outdir']+'/BAM/{sample}.bam',
        bai=config['outdir']+'/BAM/{sample}.bam.bai',
        genome=config['genome_dir']+'/hg38.fa',
        vcf=config['outdir']+'/marginPhase/{sample}.vcf'
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
        '../envs/tabix.yml'
    shell:
        'bgzip {input}'

rule final_vcf_idx:
    input:
        config['outdir']+'/VCF/{sample}.vcf.gz'
    output:
        config['outdir']+'/VCF/{sample}.vcf.gz.tbi'
    conda:
        '../envs/tabix.yml'
    shell:
        'tabix {input}'
