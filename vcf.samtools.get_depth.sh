#!/bin/bash

bname=`basename $1 .vcf.gz`

#DP for sample is alt depth 
#-GF AD is useless- not set

if [ -e $1 ]
then
    gatk-launch VariantsToTable \
     -T VariantsToTable \
     -R $2 \
     -V $1 \
     -F CHROM -F POS -F REF -F ALT  -F DP -GF DP \
     -O $bname.table
     #--allowMissingData
fi