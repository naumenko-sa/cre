#!/bin/bash

bname=`basename $1 .subset.vcf.gz`

if [ -e $1 ]
then
    gatk-launch VariantsToTable \
    -R $2 \
    -V $1 \
    -F CHROM -F POS -F REF -F ALT -F DP -GF DP -GF AD -GF GT \
    -O $bname.table
    #--allowMissingData -not used since gatk3.8
fi
