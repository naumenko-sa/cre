#!/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=30g,mem=30g

#usage: qsub ~/cre/bcl2fastq.sh -v run=<run_folder>,sample_sheet=SampleSheet.csv

#run: base run-folder where RunInfo.xml , SampleSheet.csv are present along with Data/Intensities/BaseCalls/L<lane folders>
#FASTQ output will be placed in $run/Data/Intensities/BaseCalls/
if [ -z $run ]; then
	run=$1;
fi;

#pass this if the sample_sheet name is not under run folder/named differently
if [[ -z $sample_sheet || -n $2 ]]; then
	sample_sheet=$2;
fi;

arg=" -R ${run} ";
if [ ! -z $sample_sheet ]; then 
	arg=${arg}" --sample-sheet ${sample_sheet}";
fi;



#2.19 needs 20G of RAM for 10 threads
# sometimes 20G is not enough - I've got an error
# bcl2fastq is designed up to 32G
module load bcl2fastq

#I use no-lane splitting, and the result fastq files could not be uploaded to basespace
echo "CMD: bcl2fastq $arg -p  10 --no-lane-splitting"
bcl2fastq $arg -p  10 --no-lane-splitting 

