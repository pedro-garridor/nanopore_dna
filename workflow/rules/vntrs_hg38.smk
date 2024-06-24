'''
Nanopore DNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

# NOTE: this rule only will run if user does not change default config['vntr_bed'] variable
rule get_vntr_hg38:
    output:
        vntr=temp('human_GRCh38_no_alt_analysis_set.trf.bed')
    retries: 3
    shell:
        "curl https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed > {output}"
