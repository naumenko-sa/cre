#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=15g,mem=15g

# script to cleanup the report directory after a re-run

# parameters: family = familyid 

# input: familyid 
# output: renamed re-run folder with current date

curr_date=$(date +"%Y-%m-%d")

if [ -z "$family" ]
then
	family=$1
fi

if [ ! -d $family ]
then
	echo $family " folder not present. Exiting."
	exit
fi

mv $family "${family}_${curr_date}"

# optionally, can also move the files from within subfolder to top level and/or remove logs
