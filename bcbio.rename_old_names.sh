#!/bin/bash

#renames old runs of bcbio to annotated-decomposed.vcf.gz form

if [ -z $family ];
then
    family=$1
fi

suffix="annotated-decomposed.vcf.gz"

mv ${family}-ensemble.vcf.gz ${family}-ensemble-${suffix}
mv ${family}-ensemble.vcf.gz.tbi ${family}-ensemble-${suffix}.tbi

for tool in {platypus,freebayes,gatk-haplotype};
do
    mv ${family}-${tool}.decomposed.vcf.gz ${family}-${tool}-${suffix}
    mv ${family}-${tool}.decomposed.vcf.gz.tbi ${family}-${tool}-${suffix}.tbi
done
