#!/bin/bash

#fixes sample names to family_sample

#$1 = family_id

cat samples.txt | awk -v fam=$1 '{print fam"_"$1}' > samples.txt.fixed
rm samples.txt
mv samples.txt.fixed samples.txt

for f in *ready.bam.bai;do mv $f `echo $f | sed s/"-ready"//`;done;
for f in *ready.bam;do mv $f `echo $f | sed s/"-ready"//`;done;
for f in *.bam;do mv $f `echo ${1}_${f}`;done;
for f in *.bam.bai;do mv $f `echo ${1}_${f}`;done;


for f in *.vcf.gz;do bcftools reheader -s samples.txt $f  > $f.reheader;done;
rm *.vcf.gz
for f in *.reheader;do mv $f `echo $f | sed s/.reheader//`;done;
for f in *.vcf.gz; do tabix $f;done;

vcf.split_multi.sh $1.vcf.gz

# rerun cre.gemini_load.sh
# rerun cre.sh
