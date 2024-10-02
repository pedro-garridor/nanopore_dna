'''
Nanopore DNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

configfile: "config/config.yaml"

rule hapdiff:
    input:
        genome=config['genome_fa'],
        vntr=config['vntr_bed'],
        hap=config['outdir']+'/HapDup/{sample}'
    output:
        temp(directory(config['outdir']+'/HapDiff/{sample}'))
    threads:
        workflow.cores/4
    singularity:
        'docker://mkolmogo/hapdiff:0.9'
    shell:
        '''
        hapdiff.py \
            --reference {input.genome} \
            --pat {input.hap}/hapdup_dual_1.fasta \
            --mat {input.hap}/hapdup_dual_2.fasta \
            --tandem-repeats {input.vntr} \
            --out-dir {output} \
            -t {workflow.cores}
        '''
