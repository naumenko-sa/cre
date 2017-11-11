#!/bin/bash

bname=`basename $1 .vcf.gz`

if [ -e $1 ]
then
    gatk -R ~/work/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -T VariantsToTable -V $1 \
     -F CHROM -F POS -F REF -F ALT -F DP -GF DP -GF AD -GF GT -o $bname.table 
    #--allowMissingData -not used since gatk3.8
fi