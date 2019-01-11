#!/bin/bash
#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

#variant bam: https://github.com/walaj/VariantBam
#needs to install xz and newer gcc

module load gcc/5.2.0

XZ_LIB=/hpf/largeprojects/ccmbio/naumenko/tools/xz-install/lib/
VARIANT_BAM_PATH=/hpf/largeprojects/ccmbio/naumenko/tools/VariantBam

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$XZ_LIB

bname=`basename $bam .bam`
#removing reads with grep is necessary for manta: https://github.com/walaj/VariantBam/issues/16
$VARIANT_BAM_PATH/variant $bam -L /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37d5/seq/decoy.bed -v | grep -v hs37d5 | grep -v NC_007605 | samtools view - -hb > $bname.no_decoy_reads.bam



