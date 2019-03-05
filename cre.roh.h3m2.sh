#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

#ROH analysis with 
#https://sourceforge.net/projects/h3m2/files/
#https://www.ncbi.nlm.nih.gov/pubmed/24966365

export H3M2PATH=/hpf/largeprojects/ccmbio/naumenko/tools/H3M2Tool

$H3M2PATH/H3M2BamParsing.sh $H3M2PATH . $bam . LOH_analysis `basename $bam .bam` /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa $H3M2PATH/SNP1000GP.HGb37_Exome.bed

$H3M2PATH/H3M2Analyze.sh $H3M2PATH . LOH_analysis `basename $bam .bam` $H3M2PATH/SNP1000GP.HGb37_Exome.bed 100000 0.1 0.1 5

#DNorm -> parameter of the H3M2 (it can be 1000, 10000 and 100000. We suggest to use 100000)
#P1 -> parameter of the H3M2 (set to 0.1)
#P2 -> parameter of the H3M2 (set to 0.1)
#Factor -> parameter of the H3M2 (set to 5)

