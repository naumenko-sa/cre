#!/bin/bash
# prepares full chomosome test set for gnomad to test ggd recipe
mkdir txtmp
cd txtmp
prefix=gnomad.exomes.r2.1.sites.grch38.chr
for chrom in $(seq 1 22;echo X Y)
do
  curl -r 0-1000000 -O http://ftp.ensemblorg.ebi.ac.uk/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/r2.1/exomes/${prefix}${chrom}_noVEP.vcf.gz
  gunzip -c ${prefix}${chrom}_noVEP.vcf.gz | head -n 1400 > ${prefix}${chrom}_noVEP.vcf
  bgzip -f ${prefix}${chrom}_noVEP.vcf
  tabix ${prefix}${chrom}_noVEP.vcf.gz
done
cd ..
