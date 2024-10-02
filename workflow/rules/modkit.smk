'''
Nanopore DNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

configfile: "config/config.yaml"

rule modkit:
    input:
        hap_bam=config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.haplotagged.bam',
        hap_bai=config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.haplotagged.bam.bai',
        reference=config['genome_fa']
    output:
        bed_methyl=protected(config['outdir']+'/modkit/{sample}.bed')
    conda:
        '../envs/modkit.yml'
    shell:
        '''
        modkit pileup {input.hap_bam} {output} \
            --cpg \
            --ref {input.reference}
        '''
