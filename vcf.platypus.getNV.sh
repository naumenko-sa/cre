#!/bin/bash

bname=`basename $1 | .vcf.gz`

if [ -e $1 ]
then
    gatk -R ~/work/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -T VariantsToTable -V $1 \
	 -F CHROM -F POS -F REF -F ALT -F TC -GF NR -GF NV -o $bname.table 
     #--allowMissingData 
fi
