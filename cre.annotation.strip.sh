#!/bin/bash

# also vt rminfo but it requires exact case matching

bname=`basename $1 .vcf.gz`

bcftools annotate -x INFO/CSQ,INFO/ANN,INFO/af_adj_exac_afr,INFO/af_adj_exac_amr,INFO/af_adj_exac_eas,INFO/af_adj_exac_fin,INFO/af_adj_exac_nfe,INFO/af_adj_exac_oth,INFO/af_adj_exac_sas,\
INFO/af_exac_all,INFO/ac_exac_all,INFO/an_exac_all,INFO/ac_adj_exac_afr,INFO/an_adj_exac_afr,INFO/ac_adj_exac_amr,INFO/an_adj_exac_amr,INFO/ac_adj_exac_eas,INFO/an_adj_exac_eas,\
INFO/ac_adj_exac_fin,INFO/an_adj_exac_fin,INFO/ac_adj_exac_nfe,INFO/an_adj_exac_nfe,INFO/ac_adj_exac_oth,INFO/an_adj_exac_oth,INFO/ac_adj_exac_sas,INFO/an_adj_exac_sas,\
INFO/common_pathogenic,INFO/max_aaf_all,INFO/num_exac_Het,INFO/num_exac_Hom,INFO/rs_ids,INFO/af_1kg_amr,INFO/af_1kg_eas,INFO/af_1kg_sas,INFO/af_1kg_afr,INFO/af_1kg_eur,INFO/af_1kg_all,\
INFO/fitcons,INFO/encode_consensus_gm12878,INFO/encode_consensus_h1hesc,INFO/encode_consensus_helas3,INFO/encode_consensus_hepg2,INFO/encode_consensus_huvec,\
INFO/encode_consensus_k562,INFO/dgv,INFO/hapmap1,INFO/hapmap2,INFO/gnomAD_AF,INFO/gnomAD_AFR_AF,INFO/gnomAD_AMR_AF,INFO/gnomAD_ASJ_AF,INFO/gnomAD_EAS_AF,INFO/gnomAD_FIN_AF,\
INFO/gnomAD_NFE_AF,INFO/gnomAD_OTH_AF,INFO/gnomad_ac_es,INFO/gnomad_hom_es,INFO/gnomad_af_es,INFO/gnomad_an_es,INFO/gnomad_ac_gs,INFO/gnomad_hom_gs,INFO/gnomad_af_gs,\
INFO/gnomad_an_gs,INFO/gnomad_ac,INFO/gnomad_af_popmax,INFO/gnomad_an,INFO/gnomad_af,INFO/clinvar_pathogenic,INFO/clinvar_sig,INFO/CADD_phred,INFO/phyloP20way_mammalian,\
INFO/phastCons20way_mammalian,INFO/Vest3_score,INFO/Revel_score,INFO/Gerp_score,INFO/vcfanno_gnomad_ac_es,INFO/vcfanno_gnomad_hom_es,INFO/vcfanno_gnomad_af_es,\
INFO/vcfanno_gnomad_an_es,INFO/vcfanno_gnomad_ac_gs,INFO/vcfanno_gnomad_hom_gs,INFO/vcfanno_gnomad_af_gs,INFO/vcfanno_gnomad_an_gs,INFO/vcfanno_gnomad_ac,\
INFO/vcfanno_gnomad_af_popmax,INFO/vcfanno_gnomad_an,INFO/vcfanno_gnomad_af,INFO/gnomad_gc,INFO/gnomad_gc_female,INFO/gnomad_gc_male\
	     $1 -o $bname.no_anno.vcf.gz

tabix $bname.no_anno.vcf.gz
