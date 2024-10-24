'''
Nanopore DNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

configfile: "config/config.yaml"

rule pepper:
    input:
        bam=config['outdir']+'/BAM/{sample}.bam',
        bai=config['outdir']+'/BAM/{sample}.bam.bai',
        genome=config['ref']
    output:
        intermediate_files=temp(directory(config['outdir']+'/PMDV/{sample}/intermediate_files')),
        logs=protected(directory(config['outdir']+'/PMDV/{sample}/logs')),
        chunks=temp(config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.chunks.csv'),
        hap_bam=protected(config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.haplotagged.bam'),
        vcf=protected(config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.phased.vcf.gz'),
        tbi=protected(config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.phased.vcf.gz.tbi'),
        pahset_bed=protected(config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.phaseset.bed'),
        vcf_unphased=temp(config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.vcf.gz'),
        tbi_unphased=temp(config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.vcf.gz.tbi'),
        report=protected(config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.visual_report.html')
    params:
        outdir = config['outdir']
    resources:
        nvidia_gpu=1
    threads:
        workflow.cores/2
    singularity:
        /opt/pepper_deepvariant_r0-8-gpu/ # 'docker://kishwars/pepper_deepvariant:r0.8-gpu'
    shell:
        '''
        run_pepper_margin_deepvariant \
            call_variant \
            -b {input.bam} \
            -f {input.genome} \
            -o {params.outdir}/PMDV/{wildcards.sample} \
            -t {threads} \
            -s {wildcards.sample} \
            --ont_r10_q20 \
            --phased_output \
            -p {wildcards.sample}_PMDV_FINAL
        '''

rule pepper_samtools_index:
    input:
        config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.haplotagged.bam',
    output:
        protected(config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.haplotagged.bam.bai')
    conda:
        '../envs/samtools.yml'
    shell:
        'samtools index {input}'
