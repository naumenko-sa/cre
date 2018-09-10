#!/bin/bash

bname=`basename $1 .subset.vcf.gz`

if [ -e $1 ] && [ -e $2 ]
then
    gatk VariantsToTable \
    -R $2 \
    -V $1 \
    -F CHROM -F POS -F REF -F ALT -F DP -GF DP -GF AO \
    -O $bname.table
     #--allowMissingData 
fi
