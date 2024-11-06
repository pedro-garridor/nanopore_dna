#!/bin/bash

VERSION='0.0.1+20241030'

echo ''
echo 'nanopore_dna: DNA-Seq nanopore pipeline'
echo 'Copyright (C) 2024, Pedro Garrido Rodríguez'
echo ''
echo 'This program is distributed in the hope that it will be useful,'
echo 'but WITHOUT ANY WARRANTY; without even the implied warranty of'
echo 'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.'
echo ''
echo "Version $VERSION"
echo ''

if [ "$#" -eq 0 ]; then
    echo "No arguments supplied."
    echo "Use 'bash nanopore_dna.sh' to get help."
    exit 2
fi

while [ $# -gt 0 ]; do
    case $1 in 
        -i|--input)
            INPUT="$2"
            shift
            shift
            ;;
        -o|--output)
            OUTDIR="$2"
            shift
            shift
            ;;
        -r|--reference)
            REF="$2"
            shift
            shift
            ;;
        -g|--hg38)
            HG38="$2"
            shift
            shift
            ;;
        -n|--dryrun)
            DRYRUN=1
            shift
            ;;
        -t|--threads)
            THREADS="$2"
            shift
            shift
            ;;
        -v|--version)
            exit
            shift
            ;;
        -h|--help)
            echo 'Usage: bash nanopore_dna.sh [arguments]'
            echo ' '
            echo 'Required:'
            echo '  -i, --input FILE                        raw data folder'
            echo '  -o, --output FILE                       output folder'
            echo '  -r , reference FILE                     reference genome'
            echo ' '
            echo 'Optional arguments:'
            echo '  -g, --hg38 [no]                         reference is hg38? [yes|no]'
            echo '  -t [1]                                  threads'
            echo '  -n                                      perform a dry run'
            echo ' '
            echo '[i] Please note you will need a CUDA-compatible GPU for this software.'
            exit
            ;;
        -*|--*)
            echo "Unknown option $1"
            exit 1
            ;;
        *)
            echo "Unknown argument $1"
            exit 1
            ;;
    esac
done

if [ -z $INPUT ] || [ -z $OUTDIR ] || [ -z $REF ]; then
    echo "Arguments -i, -o and -r are required."
    echo "Use 'bash nanopore_dna.sh -h' to get help."
    exit 2
fi

if [ $DRYRUN ]; then
    snakemake -n \
        -kj $THREADS \
        --config \
        input=$INPUT \
        outdir=$OUTDIR \
        ref=$REF \
        hg38=$HG38 \
        --use-conda \
        --use-singularity \
        --singularity-args "-B $INPUT,$OUTDIR,$(dirname $REF)"
else
    snakemake \
        -kj $THREADS \
        --config \
        input=$INPUT \
        outdir=$OUTDIR \
        ref=$REF \
        hg38=$HG38 \
        --use-conda \
        --use-singularity \
        --singularity-args "-B $INPUT,$OUTDIR,$(dirname $REF)"
fi
