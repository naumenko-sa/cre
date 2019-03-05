#!/bin/bash
# runs of HOM MAF 5% variant in proband (CH0620) affected sib (CH0621) and in both
# gt_types=0 or 3 = HOM, 1 = HET

sample=$1
echo "chrom,pos,ref,alt,impact,qual,dp, gene, maf, gts."$sample", gt_types."$sample",stretch"

gemini query -q "select chrom,start+1 as pos, ref, alt,impact,qual,dp, gene, gnomad_af_popmax as maf, gts."$sample",gt_types."$sample" from variants where
type='snp' and dp>20 and qual>100 and gnomad_af_popmax<0.05" --gt-filter "gt_types."$sample" != 2" $2 | sort -k1,1n -k2,2n \
| awk -F "\t" 'BEGIN{prev=1;stretch=0}{if($11 == 3){genotype=0}else{genotype=$11};if(genotype==prev && genotype==0){stretch=stretch+1}else{stretch=0};prev=genotype;print $0"\t"stretch}' | sed s/"\t"/","/g
