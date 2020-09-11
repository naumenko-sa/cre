#!/bin/bash
# generate both the regular and synonymous reports for a family
# pass in the family id as the first positional parameter, use -email to recieve an email when both jobs are done
# usage: generate_reports.sh <family id> [optional -email] 

family=$1

if [ "$2" == "-email" ]; then
  email_flag="-m e"
else
	email_flag=""
fi

cd $family
family_vcf="${family}-ensemble-annotated-decomposed.vcf.gz"
if [ -f $family_vcf ]
then
	vcf2cre_job="$(qsub ~/cre/cre.vcf2cre.sh -v original_vcf="${family}-ensemble-annotated-decomposed.vcf.gz",project=${family})"
else
	echo "${family_vcf} Not Present, Aborting."
	cd ..
	exit
fi

standard_job="$(qsub ~/cre/cre.sh -W depend=afterany:"${vcf2cre_job}" -v family=${family})"
echo "Standard Report Job ID: ${standard_job}"

synonymous_job="$(qsub ~/cre/cre.sh -W depend=afterany:"${standard_job}" -v family=${family},type=wes.synonymous ${email_flag})"
echo "Synonymous Report Job ID: ${synonymous_job}"

echo "The Rerun subfolder will be renamed by the current date after the reports are created"
cleanup_job="$(qsub ~/cre/rename_rerun.sh -W depend=afterok:"${synonymous_job}" -v family=${family})"

cd ..
