#!/bin/bash

bname=`basename $1 .vcf.gz`
bcftools annotate -x INFO/CSQ,INFO/af_adj_exac_afr,INFO/af_adj_exac_amr,INFO/af_adj_exac_eas,INFO/af_adj_exac_fin,INFO/af_adj_exac_nfe,INFO/af_adj_exac_oth,INFO/af_adj_exac_sas,INFO/af_exac_all,INFO/common_pathogenic,INFO/max_aaf_all,INFO/num_exac_Het,INFO/num_exac_Hom \
    $1 \
    -o $bname.no_anno.vcf.gz -Oz
tabix $bname.no_anno.vcf.gz

