#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

# reference should be decompressed and indexed: using one from bcbio
# usually crams come with hg19 reference (chr), not GRCh37

# parameters: cram = file.cram
#	      sample = familyid_sampleid

# input: file.cram
# output: sample_1.fq.gz sample_2.fq.gz

module load java

cramtools fastq -Xmx10g -F $sample --skip-md5-check \
		-z \
		-I $cram \
		-R /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa

mv $sample_1.fastq.gz $sample_1.fq.gz
mv $sample_2.fastq.gz $sample_2.fq.gz
