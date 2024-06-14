'''
Nanopore DNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

rule get_vntr_hg38:
    output:
        vntr=temp('rmsk.txt.gz')
    retries: 3
    shell:
        "curl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz > {output}"

rule bed_vntr_hg38:
    input:
        vntr='rmsk.txt.gz'
    output:
        bed=temp('vntrs_hg38.bed')
    conda:
        '../envs/bedtools.yml'
    shell:
        """
        zcat {input} |\
        awk '{print $6 "\t" $7 "\t" $8}' |\
        tail -n +2 |\
        bedtools sort > {output}
        """
