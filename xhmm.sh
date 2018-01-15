#!/bin/bash

#PBS -l walltime=240:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

#http://atgu.mgh.harvard.edu/xhmm/tutorial.shtml#params_file

#bam.list
#EXOME.interval_list

reference=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
xhmm_home=/hpf/largeprojects/ccmbio/naumenko/tools/xhmm/statgen-xhmm-cc14e528d909

gatk -Xmx3072m \
-T DepthOfCoverage -I bam.list -L EXOME.interval_list \
-R $reference \
-dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable \
--minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 \
--includeRefNSites \
--countType COUNT_FRAGMENTS \
-o DATA.RD.txt

${xhmm_home}/xhmm --mergeGATKdepths -o DATA.RD.txt --GATKdepths DATA.RD.txt.sample_interval_summary

gatk -Xmx3072m \
-T GCContentByInterval -L EXOME.interval_list \
-R $reference \
-o DATA.locus_GC.txt

cat DATA.locus_GC.txt | awk '{if ($2 < 0.1 || $2 > 0.9) print $1}' > extreme_gc_targets.txt

${xhmm_home}/sources/scripts/interval_list_to_pseq_reg EXOME.interval_list > EXOME.targets.reg
module load plinkseq

pseq . loc-load --locdb EXOME.targets.LOCDB --file EXOME.targets.reg --group targets --out EXOME.targets.LOCDB.loc-load --noweb

pseq . loc-stats --locdb EXOME.targets.LOCDB --group targets --seqdb seqdb --noweb | awk '{if (NR > 1) print $_}' | sort -k1 -g | awk '{print $10}' | paste EXOME.interval_list - | awk '{print $1"\t"$2}' > DATA.locus_complexity.txt

cat DATA.locus_complexity.txt | awk '{if ($2 > 0.25) print $1}' > low_complexity_targets.txt

${xhmm_home}/xhmm --matrix -r ./DATA.RD.txt --centerData --centerType target \
-o ./DATA.filtered_centered.RD.txt \
--outputExcludedTargets ./DATA.filtered_centered.RD.txt.filtered_targets.txt \
--outputExcludedSamples ./DATA.filtered_centered.RD.txt.filtered_samples.txt \
--excludeTargets ./extreme_gc_targets.txt --excludeTargets ./low_complexity_targets.txt \
--minTargetSize 10 --maxTargetSize 10000 \
--minMeanTargetRD 10 --maxMeanTargetRD 500 \
--minMeanSampleRD 25 --maxMeanSampleRD 200 \
--maxSdSampleRD 150

#Runs PCA on mean-centered data:
${xhmm_home}/xhmm --PCA -r DATA.filtered_centered.RD.txt --PCAfiles DATA.RD_PCA

#Normalizes mean-centered data using PCA information:
${xhmm_home}/xhmm --normalize -r DATA.filtered_centered.RD.txt --PCAfiles DATA.RD_PCA \
--normalizeOutput ./DATA.PCA_normalized.txt \
--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7

#Filters and z-score centers (by sample) the PCA-normalized data:
${xhmm_home}/xhmm --matrix -r ./DATA.PCA_normalized.txt --centerData --centerType sample --zScoreData \
-o ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
--outputExcludedTargets ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--outputExcludedSamples ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
--maxSdTargetRD 30

#Filters original read-depth data to be the same as filtered, normalized data:
${xhmm_home}/xhmm --matrix -r ./DATA.RD.txt \
--excludeTargets ./DATA.filtered_centered.RD.txt.filtered_targets.txt \
--excludeTargets ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--excludeSamples ./DATA.filtered_centered.RD.txt.filtered_samples.txt \
--excludeSamples ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
-o ./DATA.same_filtered.RD.txt

#Discovers CNVs in normalized data:
${xhmm_home}/xhmm --discover -p ${xhmm_home}/params.txt \
-r DATA.PCA_normalized.filtered.sample_zscores.RD.txt -R ./DATA.same_filtered.RD.txt \
-c DATA.xcnv -a DATA.aux_xcnv -s ./DATA

#Genotypes discovered CNVs in all samples:
${xhmm_home}/xhmm --genotype -p ${xhmm_home}/params.txt \
-r ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt -R ./DATA.same_filtered.RD.txt \
-g ./DATA.xcnv -F $reference \
-v ./DATA.vcf
