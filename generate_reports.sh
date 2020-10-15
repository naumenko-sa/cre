#!/bin/bash
# generate both the regular and synonymous reports for a family
# pass in the family id as the first positional parameter, use -email to recieve an email when both jobs are done
# usage: generate_reports.sh <family id> <report_type> [optional -email] 

family=$1
report_type=$2
curr_date=$(date +"%Y-%m-%d")

rerun_folder="${family}_${curr_date}"

if [ "$3" == "-email" ]; then
  email_flag="-m e"
else
	email_flag=""
fi

cd $family

if [ "$report_type" = "wes" ]; then
	family_vcf="${family}-ensemble-annotated-decomposed.vcf.gz"
	eval script="~/cre/cre.vcf2cre.sh"
elif [ "$report_type" = "wgs" ]; then
	family_vcf="${family}-gatk-haplotype-annotated-decomposed.vcf.gz"
	eval script="~/crg/crg.vcf2cre.sh"
else
	echo "Please enter a report type (either wes or wgs)"
	cd ..
	exit
fi

if [ -f $family_vcf ]; then
	mkdir $rerun_folder
	cd $rerun_folder

	# could link instead of copying but not sure whether there
	# would be side effects on the linked file (#TODO: test)
	cp ../${family_vcf} .
	vcf2cre_job="$(qsub "${script}" -v original_vcf="${family_vcf}",project=${family})"
else
	echo "${family_vcf} not present, exiting."
	cd ..
	exit
fi

if [ "$report_type" = "wes" ]; then
	standard_job="$(qsub ~/cre/cre.sh -W depend=afterany:"${vcf2cre_job}" -v family=${family})"
	echo "Standard WES Report Job ID: ${standard_job}"
elif [ "$report_type" = "wgs" ]; then
	standard_job="$(qsub ~/cre/cre.sh -W depend=afterany:"${vcf2cre_job}" -v family=${family},type=wgs)"
  echo "WGS Report Job ID: ${standard_job}"
fi

echo "The re-run subfolder will be cleaned up after the reports are created"
cleanup_job="$(qsub ~/cre/cleanup_run.sh -W depend=afterok:"${standard_job}" -v family=${family})"

cd ../..
