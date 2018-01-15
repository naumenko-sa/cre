#!/bin/bash

#PBS -l walltime=24:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

#http://atgu.mgh.harvard.edu/xhmm/tutorial.shtml#params_file

#bam.list
#EXOME.interval_list

reference=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
xhmm_home=/hpf/largeprojects/ccmbio/naumenko/tools/xhmm/statgen-xhmm-cc14e528d909

gatk -Xmx3072m \
-T DepthOfCoverage -I bam.$n.list -L EXOME.interval_list \
-R $reference \
-dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable \
--minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 \
--includeRefNSites \
--countType COUNT_FRAGMENTS \
-o group${n}.DATA
