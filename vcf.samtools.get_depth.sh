#!/bin/bash

bname=`basename $1 .vcf.gz`

#DP for sample is alt depth 
#-GF AD is useless- not set

#find reference
reference=`which gatk-launch | sed s/"bin\/gatk-launch"/"bcbio\/genomes\/Hsapiens\/GRCh37\/seq\/GRCh37.fa"/`

if [ -e $1 ]
then
    gatk-launch VariantsToTable \
     -T VariantsToTable \
     -V $1 \
     -F CHROM -F POS -F REF -F ALT  -F DP -GF DP \
     -O $bname.table
     #--allowMissingData
fi