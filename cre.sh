#!/bin/bash

# cleans up after bcbio - 
# 	when running large cohort only the final folder is kept and only one ensemble gemini database: 2-3G per family
# 	keeps bam files for new samples
# prepares tables for report generation
# generates report

# parameters:
# family = [family_id] (=folder_name,main result file should be family-ensemble.db)
# cleanup= [0|1] default = 0
# make_report=[0|1] default = 1
# rnaseq=[0|1] default = 0

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g


function f_cleanup
{

    # better to look for project-summary than hardcode the year
    # keep bam files for new samples
    cd $family
    result_dir=`find final -name project-summary.yaml | sed s/"\/project-summary.yaml"//`
    if [ -d $result_dir ];
    then
	mv $result_dir/* .
	mv final/*/*.bam .
	mv final/*/*.bai .
        rm -rf final/
	rm -rf work/
    fi

    #don't remove input files for new projects
    #rm -rf input/

    #rename bam files to match sample names
    for f in *ready.bam;do mv $f `echo $f | sed s/"-ready"//`;done;
    for f in *ready.bam.bai;do mv $f `echo $f | sed s/"-ready"//`;done;

    #make bam files read only
    for f in *.bam;do chmod 444 $f;done;

    #calculate md5 sums
    for f in *.bam;do md5sum $f > $f.md5;done;

    #validate bam files
    for f in *.bam;do	cre.bam.validate.sh $f;done;
    
    # we don't need gemini databases for particular calling algorythms
    rm ${family}-freebayes.db
    rm ${family}-gatk-haplotype.db
    rm ${family}-samtools.db
    rm ${family}-platypus.db

    cd ..
}

function f_make_report
{
    cd $family

    cre.gemini2txt.sh ${family}-ensemble.db $depth_threshold $severity_filter
    cre.gemini_variant_impacts.sh ${family}-ensemble.db $depth_threshold $severity_filter

    for f in *.vcf.gz;
    do
	tabix $f;
    done

    # report filtered vcf for import in phenotips
    # note that if there is a multiallelic SNP, with one rare allele and one frequent one, both will be reported in the VCF,
    # and just a rare one in the excel report
    cat ${family}-ensemble.db.txt | cut -f 23,24  | sed 1d | sed s/chr// > ${family}-ensemble.db.txt.positions
    bcftools view -R ${family}-ensemble.db.txt.positions -o ${family}.vcf.gz -O z ${family}-ensemble-annotated-decomposed.vcf.gz
    tabix $family.vcf.gz

    #individual vcfs for uploading to phenome central
    vcf.split_multi.sh $family.vcf.gz

    vcf.ensemble.getCALLERS.sh $family.vcf.gz

    #decompose first for the old version of bcbio!
    #gemini.decompose.sh ${family}-freebayes.vcf.gz
    vcf.freebayes.getAO.sh ${family}-freebayes-annotated-decomposed.vcf.gz

    #gemini.decompose.sh ${family}-gatk-haplotype.vcf.gz
    vcf.gatk.get_depth.sh ${family}-gatk-haplotype-annotated-decomposed.vcf.gz

    #gemini.decompose.sh ${family}-platypus.vcf.gz
    vcf.platypus.getNV.sh ${family}-platypus-annotated-decomposed.vcf.gz

    cd ..

    # using Rscript from bcbio
    Rscript ~/cre/cre.R $family
    
    cd $family
    rm $family.create_report.csv $family.merge_reports.csv
    cd ..
}

if [ -z $family ]
then
    family=$1
fi

echo $family

if [ -z "$rnaseq" ];
then
    rnaseq=0
fi

export depth_threshold=10
if [ "$rnaseq" -eq 1 ];
then
    export depth_threshold=5
    export severity_filter=ALL
fi

#no cleanup by default
if [ -z $cleanup ]
then
    cleanup=0
fi

if [ $cleanup -eq 1 ]
then
    f_cleanup
fi

if [ -z $make_report ]
then
    make_report=1
fi 

if [ $make_report -eq 1 ]
then
    f_make_report
fi
