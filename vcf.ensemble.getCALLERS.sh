#!/bin/bash

bname=`basename $1 .vcf.gz`

gatk-launch VariantsToTable \
    -R $2 -V $1 \
    -F CHROM -F POS -F REF -F ALT -F CALLERS \
    -O $bname.table

#--allowMissingData 
