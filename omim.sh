#!/bin/bash

cat genemap2.txt | grep -v "^#" | grep '(3)' | grep ENSG | awk -F "\t" '{print $11"\t"$13}' | sort -k 1,1 > omim.tmp

#some genes may have two entries

cat omim.tmp  | awk -F "\t" 'BEGIN{prev_gene="Ensembl_gene_id\tOmim_gene_description";buf=""}{if(prev_gene != $1){print prev_gene"\t"buf;buf=$2;prev_gene=$1}else{buf=buf","$2;}}END{print prev_gene","buf'} > omim.txt

rm omim.tmp