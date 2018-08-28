#!/bin/bash

bcbio_setup_genome.py -f $1 -g /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/rnaseq/ref-transcripts.gtf -n Hsapiens -b GRCh37d5 -i bwa star rtg

