#!/bin/bash
# use case:  when somebody wants to see all variants in a gene
# usage: 
# cre.gemini.get_variants4gene.sh [db] [sample] [gene]
# output:
# gene.csv

db=$1
sample=$2
gene=$3

query="select chrom,start+1 as pos,ref,alt,type,sub_type,gene,exon,codon_change,impact,gts."$sample",dp,rs_ids,gnomad_af_popmax,clinvar_sig,cadd_phred,hgvsc from variants where gene='"$gene"'"
echo $query
gemini query --header -q "$query" $db | sed s/"\t"/","/g > $gene.csv

