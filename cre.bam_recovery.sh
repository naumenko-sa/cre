#!/bin/bash

#PBS -l walltime=5:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

#http://genome.sph.umich.edu/wiki/BamUtil:_convert#BAM_File_Recovery

bam convert --recover --in $bam --out `echo $bam | sed s/bam/recovered.bam/`
