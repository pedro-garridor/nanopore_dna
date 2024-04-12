'''
Nanopore DNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

rule pepper:
    input:
        bam=config['outdir']+'/BAM/{sample}.bam',
        bai=config['outdir']+'/BAM/{sample}.bam.bai',
        genome=config['genome_dir']+'/hg38.fa'
    output:
        temp(directory(config['outdir']+'/PMDV/{sample}'))
    resources:
        nvidia_gpu=1
    threads:
        workflow.cores/4
    singularity:
        'docker://kishwars/pepper_deepvariant:r0.8-gpu'
    shell:
        '''
        run_pepper_margin_deepvariant \
            call_variant \
            -b {input.bam} \
            -f {input.genome} \
            -o {output} \
            -t {threads} \
            -s {wildcards.sample} \
            --ont_r10_q20 \
            --phased_output \
            -p {wildcards.sample}_PMDV_FINAL
        '''

