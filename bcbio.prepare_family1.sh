#!/bin/bash

#prepare family for bcbio run when input files are family_sample.bam or family_sample_1/2.fq.gz

family=$1

cd $family

cp ~/cre/bcbio.sample_sheet_header.csv $family.csv

cd input

ls | sed s/.bam// | sed s/.bai// | sed s/"_1.fq.gz"// | sed s/"_2.fq.gz"// | sort | uniq > ../samples.txt


cd ..

while read sample
do
    echo $sample","$sample","$family",,," >> $family.csv
done < samples.txt

bcbio_nextgen.py -w template ~/cre/bcbio.templates.exome.yaml $family.csv input/*

mkdir config
mkdir work

cp $family/config/$family.yaml config

cd ..


