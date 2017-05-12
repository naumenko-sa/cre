#!/bin/bash

# creates gemini database with vep refseq annotations, 
# dumps text file and variant_impacts for rare harmful variants
# it is easier just to get variant super index - VEP refseq annotation from VCF file
# finally no databaes for refseq, but possible to create in future if needed

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

if [ -z $family ];
then
    family=$1
fi

# 1.remove VEP ensembl annotation, the file is already decomposed
vt rminfo ${family}-ensemble-annotated-decomposed.vcf.gz -t CSQ -o ${family}.no_vep.vcf.gz

# annotate with refseq - lets annotate only those variants we really need - much faster
#gemini.vep.refseq.sh ${family}.no_vep.vcf.gz

# 2.extract positions of the variants in the final report
cat ${family}-ensemble.db.txt  | awk -F "\t" '{print $23"\t"$24}' | sed 1d | sed s/chr// > $family.report.tab

# 3.extract variants
bcftools view -o $family.no_vep.rare_only.to_sort.vcf -R $family.report.tab $family.no_vep.vcf.gz

# 4.sort vcf
picard SortVcf I=$family.no_vep.rare_only.to_sort.vcf O=$family.no_vep.rare_only.vcf
bgzip $family.no_vep.rare_only.vcf
tabix $family.no_vep.rare_only.vcf.gz

# 5.annotate with refseq
gemini.vep.refseq.sh $family.no_vep.rare_only.vcf.gz

# 6.convert vcf to table
gatk -R ~/work/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -T VariantsToTable -V $family.no_vep.rare_only.vepeffects_refseq.vcf.gz -F CHROM -F POS -F REF -F ALT -F CSQ -o $family.refseq.table -U ALLOW_SEQ_DICT_INCOMPATIBILITY

# 7.parse vcf and VEP CSQ field
gemini.vep.parse.pl $family.refseq.table > $family.refseq.txt

rm $family.report.tab $family.tosort.vcf $family.refseq.table

#gemini.vep2gemini.sh ${family}.no_vep.decomposed.vepeffects_refseq.vcf.gz

#gemini.gemini2txt.sh ${family}.no_vep.decomposed.vepeffects_refseq.db

#gemini.variant_impacts.sh ${family}.no_vep.decomposed.vepeffects_refseq.db
