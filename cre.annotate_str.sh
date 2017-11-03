#!/bin/bash

#use cre.annotate_str.sh vcf_from_lobstr.vcf

bname=`basename $1 .vcf`

cat $1 | grep "^#" > $bname.sorted.vcf
#sometimes lobSTR gives you unsorted vcf for some reason
cat $1 | grep -v "^#" | sort -t $'\t' -k1,1 -k2,2n >> $bname.sorted.vcf

bgzip $bname.sorted.vcf 
tabix $bname.sorted.vcf.gz 

#there may be no gene
vcfanno -base-path ~/cre/ ~/cre/cre.annotate_str.toml $bname.sorted.vcf.gz > $bname.annotated.vcf
rm $bname.insertions.txt
echo -e "CHR\tPOS\tREF\tALT\tGENE\tGENOTYPE\tINSERTION_LENGTH_BP" > $bname.insertions.txt
cat $bname.annotated.vcf | grep -v "^#" | sed s/"gene="/"\t"/g | grep "0/1:" | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$9"\t"$NF}' | sed s/":.*"// | \
    awk '{if (length($3)<length($4)) print $0"\t"length($4)-length($3);}' | sort -k7,7 -rn >> $bname.insertions.txt
