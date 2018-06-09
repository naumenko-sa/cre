#!/bin/bash

bname=`basename $1 .vcf.gz`

#find reference
reference=`which gatk-launch | sed s/"bin\/gatk-launch"/"bcbio\/genomes\/Hsapiens\/GRCh37\/seq\/GRCh37.fa"/`

if [ -e $1 ]
then
    gatk-launch VariantsToTable \
    -R $reference \
    -V $1 \
    -F CHROM -F POS -F REF -F ALT -F TC -GF NR -GF NV \
    -O $bname.table 
     #--allowMissingData 
fi
