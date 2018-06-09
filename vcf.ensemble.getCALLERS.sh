#!/bin/bash

bname=`basename $1 .vcf.gz`

#find reference
reference=`which gatk-launch | sed s/"bin\/gatk-launch"/"bcbio\/genomes\/Hsapiens\/GRCh37\/seq\/GRCh37.fa"/`

gatk-launch VariantsToTable \
    -R $reference -V $1 \
     -F CHROM -F POS -F REF -F ALT -F CALLERS -o $bname.table

#--allowMissingData 
