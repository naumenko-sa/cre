#!/bin/bash

# input is a text file of newline separated family ids
family_ids=$1
echo $family_ids

if [ "$2" = "-dryrun" ]; then
	dryrun=true
else
  dryrun=false
fi

while read family;
do
	# if the family already isn't in the current folder
	if [ ! -d $family ]
	then
		echo "Searching for " $family " in Results"
		family_folder=$(ls -d /hpf/largeprojects/ccm_dccforge/dccforge/results/*/$family 2>/dev/null;)
		if [ ! -z $family_folder ]
		then
			ln -s $family_folder .
		else
			echo "Family " $family " Not Found"
		fi		
	else
		echo "Folder for " $family " Already Exists"
	fi

	# if it was successful or already existed, list all the current reports
	if [ -d $family ]
	then
		echo "Reports:"
		ls $family/*20*.csv;
	fi

	# submit a job to regenerate (uncomment to run)
  if [ "$dryrun" = false ]
	then
	  ~/cre/generate_reports.sh $family # generate both report types
  else
    echo "dryrun set, will not generate reports for ${family}"
  fi

done < ${family_ids}
