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

standard_job="$(qsub ~/cre/cre.sh -v family=${family})"
echo "Standard Report Job ID: ${standard_job}"

synonymous_job="$(qsub ~/cre/cre.sh -W depend=afterany:"${standard_job}" -v family=${family},type=wes.synonymous ${email_flag})"
echo "Synonymous Report Job ID: ${synonymous_job}"
