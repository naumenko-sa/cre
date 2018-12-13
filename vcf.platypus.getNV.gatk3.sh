#!/bin/bash

# if vcf does not have contig lengths in the header, gatk4 fails

bname=`basename $1 .subset.vcf.gz`

if [ -e $1 ]
then
    gatk3 -Xmx10G -T VariantsToTable \
    -R $2 \
    -V $1 \
    -F CHROM -F POS -F REF -F ALT -F TC -GF NR -GF NV \
    -o $bname.table 
    #--allowMissingData 
fi
