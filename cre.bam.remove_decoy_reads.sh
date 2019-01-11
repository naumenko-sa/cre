#!/bin/bash
#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

module load gcc/5.2.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/hpf/largeprojects/ccmbio/naumenko/tools/xz-install/lib/

/hpf/largeprojects/ccmbio/naumenko/tools/VariantBam/variant 159_CH0315-ready.bam -L /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37d5/seq/decoy.bed -v | grep -v hs37d5 | grep -v NC_007605 | samtools view - -hb > 159_CH0315.bam
