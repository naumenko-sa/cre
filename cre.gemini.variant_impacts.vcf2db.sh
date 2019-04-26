#!/bin/bash
#   exports variant_impacts from gemini.db database to gemini.db.variant_impacts.txt file
#   database schema: https://gemini.readthedocs.io/en/latest/content/database_schema.html#the-variants-table
#   by default bcbio writes PASS only variants to the database

#PBS -l walltime=1:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

if [ -z $file ]
then
    file=$1
fi

depth_threshold=$2

severity_threshold=$3

max_af=$4

if [[ "$severity_threshold" == 'ALL' ]]
then
#used for RNA-seq = 20k variants in the report
    severity_filter=""
#use for WES = 1k variants in the report
elif [[ "$severity_threshold" == "wes.synonymous" ]]
then
    severity_filter="(v.is_coding=1 or v.is_splicing=1) and "
else
    severity_filter="v.impact_severity<>'LOW' and "
fi
        

sQuery="select \
	i.variant_id,\
	i.gene,\
	i.transcript,\
	i.is_exonic,\
	i.is_coding,\
	i.exon,\
	i.codon_change,\
	i.aa_change,\
	i.aa_length,\
	i.biotype,\
	i.impact,\
	i.impact_so,\
	i.impact_severity,\
	i.polyphen_pred,\
	i.polyphen_score,\
	i.sift_pred,\
	i.sift_score,\
	i.ccds,\
	i.hgvsc,\
	i.hgvsp,\
	v.source"

#old runs before Oct2017 does not have maxentscanfields in the database
if gemini db_info $1 | grep -q "maxentscan";
then 
    sQuery=$sQuery",\
	i.maxentscan_alt,\
	i.maxentscan_diff,\
	i.maxentscan_ref,\
	i.spliceregion"
fi

sQuery=$sQuery" from variants v,variant_impacts i where "$severity_filter"v.gnomad_af_popmax <= "$max_af" and \
		v.variant_id=i.variant_id and \
		(v.dp>="$depth_threshold" or v.dp='' or v.dp is null)"

#echo $sQuery


s_gt_filter=''
if [ -n "$denovo" ] && [ "$denovo" == 1 ]
then
    proband=`gemini query -q "select name from samples where phenotype=2" $file`
    mom=`gemini query -q "select name from samples where phenotype=1 and sex=2" $file`
    dad=`gemini query -q "select name from samples where phenotype=1 and sex=1" $file`
    
    s_gt_filter="gt_types."$proband" == 1 and (gt_types."$dad" == 0 or gt_types."$dad" == 3) and (gt_types."$mom" == 0 or gt_types."$mom" == 3)"
#		 (gt_types."$proband" == )"
    sQuery=$sQuery" and qual>=500"
    gemini query -q "$sQuery" --gt-filter "$s_gt_filter" --header $file
else
    gemini query --header -q "$sQuery" $file
fi
