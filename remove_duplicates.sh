#!/bin/bash
####################################################################################################
#   removes duplicate reads from a bam
#   parameters:
# I = [input.bam] 
# O = [output.bam]
# M = [metrics.txt]
####################################################################################################

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

unset JAVA_HOME
export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin:$PATH
picard MarkDuplicates I=$I O=$O M=$M REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT ASSUME_SORT_ORDER=coordinate
