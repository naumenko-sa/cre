#!/bin/bash

# prepares family for bcbio run when input files are family_sample.bam or family_sample_1/2.fq.gz
family=$1

# $2 = template type, 
# default = no value = default WES
# noalign - no alignment (for rerunning), 
# fast - no realignment,recalibration, and only gatk
# validation = NA12878 validation

template_type=$2
echo $template_type

cd $family

cp ~/cre/bcbio.sample_sheet_header.csv $family.csv

cd input

shopt -s extglob
ls @(*.bam|*.gz) | sed s/.bam// | sed s/.bai// | sed s/"_1.fq.gz"// | sed s/"_2.fq.gz"// | sort | uniq > ../samples.txt

cd ..

variant_regions=""

#default template
template=~/cre/cre.bcbio.templates.wes.yaml

if [ -n "$2" ]
then
    if [ $template_type == "noalign" ]
    then
	template=~/cre/cre.bcbio.templates.wes.noalign.yaml
    elif [ $template_type == "fast" ]
    then
	echo fast
	template=~/cre/cre.bcbio.templates.wes.fast.yaml
    elif [ $template_type == "validation" ]
    then
	template=~/cre/cre.bcbio.templates.wes.validation.yaml
	variant_regions=$3
    elif [ $template_type == "gatk4" ]
    then
	template=~/cre/cre.bcbio.templates.wes.gatk4.yaml
	variant_regions=$3
    fi
fi

while read sample
do
    echo $sample","$sample","$family",,,$variant_regions" >> $family.csv
done < samples.txt


bcbio_nextgen.py -w template $template $family.csv input/*

mv $family/config .
mv $family/work .
rm $family.csv
rmdir $family

cd ..
