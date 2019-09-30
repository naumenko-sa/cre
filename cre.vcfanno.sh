#!/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

if [ -z $vcf ]
then
    vcf=$1
fi

bname=`basename $vcf .vcf.gz`

prefix=$HOME/cre

vcfanno -p 5 -lua $prefix/cre.vcfanno.lua \
	     -base-path /hpf/largeprojects/ccmbio/naumenko/tools/bcbio_1.1.5/genomes/Hsapiens/GRCh37 \
	     $prefix/cre.vcfanno.conf \
	     $vcf | sed -e 's/Number=A/Number=1/g' | bgzip -c > $bname.annotated.vcf.gz

tabix $bname.annotated.vcf.gz
