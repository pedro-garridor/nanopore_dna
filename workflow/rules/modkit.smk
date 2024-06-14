rule modkit:
    input:
        hap_bam=config['outdir']+'/PMDV/{sample}/{sample}_PMDV_FINAL.haplotagged.bam',
        reference=config['genome_dir']
    output:
        bed_methyl=protected(config['outdir']+'/modkit/{sample}.bed')
    threads:
        workflow.cores/2
    conda:
        '../envs/modkit.yml'
    shell:
        '''
        modkit pileup {input.hap_bam} {output} \
            --cpg \
            --ref {input.reference} \
            -t {threads}
        '''
