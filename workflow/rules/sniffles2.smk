'''
Nanopore DNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

configfile: "config/config.yaml"

if config['hg38'] == 'Yes' or config['hg38'] == 'yes':
    rule sniffles2:
        input:
            hap_bam=config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.haplotagged.bam',
            reference=config['ref'],
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
else:
    rule sniffles2:
        input:
            hap_bam=config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.haplotagged.bam',
            reference=config['ref']
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
                --phase
            '''


rule sniffles2_gzip:
    input:
        vcf_plain=config['outdir']+'/sniffles2/{sample}.vcf'
    output:
        vcf_zipped=protected(config['outdir']+'/sniffles2/{sample}.vcf.gz')
    conda:
        '../envs/samtools.yml'
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
        '../envs/samtools.yml'
    shell:
        '''
        tabix {input}
        '''
