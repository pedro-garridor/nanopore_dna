# nanopore_dna
 Variant calling and methylation pipeline for nanopore DNA-Seq data

`nanopore_dna` is a variant calling, haplotyping and methylation detection pipeline for nanopore DNA-Seq experiments.

This implementation allows the analysis of nanopore DNA-Seq raw data (from uBAM) either on local workstations, as well as in shared servers or clusters, **without requiring admin permissions**.

Please note you will need, however, a **CUDA-compatible GPU** to call variants. A CPU-based pipeline [has been proposed](https://github.com/pedro-garridor/nanopore_dna/issues/3), although this version would require much more time to finish, and thus may not be implemented in the near future.

# Get started

To use this pipeline, you will need to install [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html), [snakemake](https://snakemake.readthedocs.io/en/stable/) and [SingularityCE](https://github.com/sylabs/singularity/releases/). To do this (without sudo permissions):

**1.** Download in your home the [latest Miniforge distribution](https://github.com/conda-forge/miniforge/releases/latest/) selecting your OS and architecture. For example, for Linux x64 it would be:
   
   ```
   wget https://github.com/conda-forge/miniforge/releases/download/<MINIFORGE_VERSION>/Mambaforge-<MINIFORGE_VERSION>-Linux-x86_64.sh
   ```
   

Changing `<MINIFORGE_VERSION>` with the latest version available.

   If you already have mamba installed in your account, go to step **2**.

**2.** Clone this repository in your system with:

    
    git clone https://github.com/pedro-garridor/nanopore_dna.git
    

**3.** Move into `nanopore_dna` and create its mamba env:

    
    cd nanopore_dna
    mamba env create -f nanopore_dna.yml

**NOTE**: this step assumes you have installed Miniforge and mamba, instead of conda. Please note if you are using conda, you may find inconsistencies when creating this env or during runtime.
    

**4.** Activate the environment:

    mamba activate nanopore_dna

**5.** Once you've done all of the above, you can run the pipeline with

    
    bash nanopore_dna.sh \
        -i <INPUT_FOLDER> \
        -o <OUTPUT_FOLDER> \
        -r <REFERENCE_GENOME> \
        -t <THREADS>
    

where:

- `<INPUT_FOLDER>`: directory where all uBAM files are located. One per sample, named with the sample ID.
- `<OUTPUT_FOLDER>`: path where you want to get the results.
- `<REFERENCE_GENOME>`: path to your genome `.fa`. `.fai` index is required to be present on the same folder.
- `<THREADS>`: threads you want the pipeline to use.

If you need more information on how to run the pipeline, you can run:

    
    bash nanopore_dna.sh -h
    

# Citations

This pipeline is part of my PhD Thesis. If you're using it, please cite it :)

> P. Garrido-Rodríguez, “Aplicación de la bioinformática en la descripción y resolución de patologías hematológicas y mecanismos biológicos relacionados,” Ph.D dissertation, Universidad de Murcia, Spain, 2025. [https://digitum.um.es/digitum/handle/10201/152160](https://digitum.um.es/digitum/handle/10201/152160)
