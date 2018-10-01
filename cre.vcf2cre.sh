#!/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

# 50g is crucial -20,30 crashes sometimes

# annotates vcf with VEP and vcfanno for cre for cre report generation
# parameters:
# original_vcf = file.vcf.gz
# project = case = family = S11 (example)
# [optional] ped = file.ped

# qsub ~/cre/cre.vcf2cre.sh -v original_vcf=file.vcf,project=412

bname=`basename $original_vcf .vcf.gz`

echo "###############################################"
echo `date` ": Removing annotation of " $original_vcf
echo "###############################################"

cre.annotation.strip.sh $original_vcf

gunzip -c $bname.no_anno.vcf.gz | grep "^#"  > $project.vcf
gunzip -c $bname.no_anno.vcf.gz | grep -v "^#" | grep PASS | grep -v possible_rnaedit | egrep -v "^GL000" >> $project.vcf
bgzip $project.vcf
tabix $project.vcf.gz

bcftools sort -Oz $project.vcf.gz > $project.sorted.vcf.gz
tabix $project.sorted.vcf.gz

cre.vt.decompose.sh $project.sorted.vcf.gz

echo "################################################"
echo `date` "Annotating with VEP..."
echo "################################################"

cre.vep.sh $project.sorted.decomposed.vcf.gz

echo "################################################"
echo `date` "Annotating with vcfanno ..."
echo "################################################"
cre.vcfanno.sh $project.sorted.decomposed.vepeffects.vcf.gz

if [ -z $ped ]
then
    echo "no ped file, generating ..."
    bcftools query -l $original_vcf > samples.txt
    > $project.ped
    for sample in `cat samples.txt`
    do
    	echo -e "1\t"$sample"\t0\t0\t0\t0\n" >> $project.ped
    done
    ped=$project.ped
fi

echo "#################################################"
echo `date` "Generating gemini database"
echo "#################################################"
vcf2db.py $project.sorted.decomposed.vepeffects.annotated.vcf.gz 	$ped	${project}-ensemble.db

mkdir $project

mv ${project}-ensemble.db $project

mv $project.sorted.decomposed.vepeffects.annotated.vcf.gz ${project}/${project}-ensemble-annotated-decomposed.vcf.gz
mv $project.sorted.decomposed.vepeffects.annotated.vcf.gz.tbi ${project}/${project}-ensemble-annotated-decomposed.vcf.gz.tbi

cd $project
ln -s ${project}-ensemble-annotated-decomposed.vcf.gz ${project}-gatk-haplotype-annotated-decomposed.vcf.gz
ln -s ${project}-ensemble-annotated-decomposed.vcf.gz.tbi ${project}-gatk-haplotype-annotated-decomposed.vcf.gz.tbi

cd ..

rm $bname.no_anno.vcf.gz
rm $bname.no_anno.vcf.gz.tbi

rm $project.vcf.gz
rm $project.vcf.gz.tbi

rm $project.sorted*.vcf.gz
rm $project.sorted*.vcf.gz.tbi

echo "#####################################################"
echo `date` " DONE"
echo "#####################################################"
