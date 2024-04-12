'''
Nanopore DNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

rule hapdiff:
    input:
        genome=config['genome_dir']+'/hg38.fa',
        hap=config['outdir']+'/HapDup/{sample}'
    output:
        temp(directory(config['outdir']+'/HapDiff/{sample}'))
    threads:
        workflow.cores/4
    singularity:
        'docker://mkolmogo/hapdiff'
    shell:
        '''
        hapdiff.py \
            --reference {input.genome} \
            --pat {input.hap}/hapdup_dual_1.fasta \
            --mat {input.hap}/hapdup_dual_2.fasta \
            --out-dir {output} \
            -t {workflow.cores}
        '''
