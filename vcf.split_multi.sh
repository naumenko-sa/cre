#!/bin/bash

# $1 = family.vcf.gz

for sample in `bcftools queru -l $1`
do 
    bcftools view -c1 -Ov -s $sample -o $sample.vcf $1
done

