#!/bin/bash

gemini query --header -q "select chrom,start+1,ref,alt,type,sub_type,gene,exon,codon_change,impact,gts.1247R_329260,depth,rs_ids,max_aaf_all,clinvar_sig,cadd_scaled,vep_hgvsc from variants where gene='NLRP3'" 1247R-ensemble.db > NLRP3.tsv
