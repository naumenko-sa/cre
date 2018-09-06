#!/bin/bash

#PBS -l walltime=20:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

bname=`basename $vcf .vcf.gz`

vt normalize -r ~/work/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa $vcf | vt uniq - -o $bname.normalized.vcf.gz
tabix $bname.normalized.vcf.gz
