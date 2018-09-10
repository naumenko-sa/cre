#!/bin/bash

bname=`basename $1 .subset.vcf.gz`

if [ -e $1 ]
then
    gatk VariantsToTable \
    -R $2 \
    -V $1 \
    -F CHROM -F POS -F REF -F ALT -F TC -GF NR -GF NV \
    -O $bname.table 
     #--allowMissingData 
fi
