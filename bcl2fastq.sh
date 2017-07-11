#!/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=30g,mem=30g

if [ -z $sample_sheet ];
then
    sample_sheet=$1
fi

#2.19 needs 20G of RAM for 10 threads
# sometimes 20G is not enough - I've got an error
# bcl2fastq is designed up to 32G
module load bcl2fastq

#I use no-lane splitting, and the result fastq files could not be uploaded to basespace
bcl2fastq --sample-sheet $sample_sheet -p  10 --no-lane-splitting
