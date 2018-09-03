#!/bin/bash

# crashes on qlogin node

#PBS -l walltime=48:00:00,nodes=1:ppn=20
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

hostname
echo $PATH
echo $PYTHONPATH

bcbio_setup_genome.py -f hs37d5.fa -g /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/rnaseq/ref-transcripts.gtf -n Hsapiens -b GRCh37d5 -i bwa star rtg -c 20

