#!/bin/bash

# if input vcf is from TCAG (HAS) it does not have DP INFO field, we need to fake it from FORMAT DP for SNVs and from DPI for indels:

bname=`basename $1 .vcf.gz`

# remove chr
gunzip -c $1 |  sed s/"ID=chrM"/"ID=MT"/ | sed s/"^chrM"/MT/ | sed s/"ID=chr"/"ID="/ | sed s/"^chr"// > $bname.nochr.vcf
bgzip $bname.nochr.vcf
tabix $bname.nochr.vcf.gz

gunzip -c $bname.nochr.vcf.gz | grep "^#" | head -n1 > $bname.dp.vcf
echo "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">" >> $bname.dp.vcf
gunzip -c $bname.nochr.vcf.gz | grep "^#" | sed 1d >> $bname.dp.vcf

#process SNVs
gunzip -c $bname.nochr.vcf.gz | grep -v "^#"  | grep PASS | grep ":DP:" | awk -F ':' '{print $0"\tDP="$9}' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11";"$8"\t"$9"\t"$10}' >> $bname.dp.vcf

#process indels
gunzip -c $bname.nochr.vcf.gz | grep -v "^#"  | grep PASS | grep ":DPI:" | awk -F ':' '{print $0"\tDP="$8}' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11";"$8"\t"$9"\t"$10}' | sed s/":DPI:"/":DP:"/ >> $bname.dp.vcf

bgzip $bname.dp.vcf
tabix $bname.dp.vcf.gz

bcftools sort -o $bname.dp.sorted.vcf.gz -Oz $bname.dp.vcf.gz
tabix $bname.dp.sorted.vcf.gz

rm $bname.nochr.vcf.gz $bname.dp.vcf.gz
