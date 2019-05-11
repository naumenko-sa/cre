#!/bin/bash
# prepares full chomosome test set for gnomad to test ggd recipe
mkdir txtmp
cd txtmp
type="genomes"
prefix=gnomad.$type.r2.1.sites.chr
# Y is absent for gnomad genomes
for chrom in $(seq 1 22;echo X Y)
do
  curl -r 0-1000000 -O http://ftp.ensemblorg.ebi.ac.uk/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad/r2.1/${type}/${prefix}${chrom}_noVEP.vcf.gz
  gunzip -c ${prefix}${chrom}_noVEP.vcf.gz | head -n 1200 > ${prefix}${chrom}_noVEP.vcf
  bgzip -f ${prefix}${chrom}_noVEP.vcf
  tabix ${prefix}${chrom}_noVEP.vcf.gz
done
cd ..
