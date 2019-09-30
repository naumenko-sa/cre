#!/bin/bash
#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

# 50g is crucial -20,30 crashes sometimes

# load annotated vcf to gemini.db 
# parameters:
# vcf = file.vcf.gz, no tabix index is needed in the dir
# project = case = family = S11 (example)

# qsub ~/cre/cre.vcf2cre.sh -v vcf=file.vcf.gz,project=412

. /hpf/largeprojects/ccmbio/naumenko/tools/bcbio_1.1.5/.test_profile

bname=`basename $vcf .vcf.gz`

bcftools query -l $vcf > samples.txt
> $project.ped
for sample in `cat samples.txt`
do
	echo -e "1\t"$sample"\t0\t0\t0\t0" >> $project.ped
done
ped=$project.ped

vcf2db.py $vcf.vcf.gz	$ped	${project}-ensemble.db

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
