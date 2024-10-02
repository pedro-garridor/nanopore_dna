'''
Nanopore DNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

configfile: "config/config.yaml"

rule sniffles2:
    input:
        hap_bam=config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.haplotagged.bam',
        reference=config['genome_fa'],
        vntr=config['vntr_bed']
    output:
        vcf=temp(config['outdir']+'/sniffles2/{sample}.vcf')
    conda:
        '../envs/sniffles2.yml'
    shell:
        '''
        sniffles \
            --input {input.hap_bam} \
            --vcf {output} \
            --reference {input.reference} \
            --tandem-repeats {input.vntr} \
            --phase
        '''

rule sniffles2_gzip:
    input:
        vcf_plain=config['outdir']+'/sniffles2/{sample}.vcf'
    output:
        vcf_zipped=protected(config['outdir']+'/sniffles2/{sample}.vcf.gz')
    conda:
        '../envs/bcftools.yml'
    shell:
        '''
        bgzip {input}
        '''

rule sniffles2_tbi:
    input:
        vcf_zipped=config['outdir']+'/sniffles2/{sample}.vcf.gz'
    output:
        tbi=protected(config['outdir']+'/sniffles2/{sample}.vcf.gz.tbi')
    conda:
        '../envs/bcftools.yml'
    shell:
        '''
        tabix {input}
        '''
