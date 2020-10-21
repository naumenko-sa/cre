#!/bin/bash

#usage: sh gemini_compare.sh <crg2-ensemble.db> <crg-ensemble.db> <pos.bed>

#$1 and $2 are ensemble dbs to extract variants
#$3 three-column(0-start) bed intervals (this could be common/uniq variants)

#assuming the same sample naming in both the dbs use this instead of gts.(*)
get_samples () {
a=($(gemini query -q "select name from samples;" $1));
for i in ${a[*]}; do gt+=($(echo "gts.${i},gt_alt_depths.${i},gt_depths.${i}")); done;
gt=$(IFS=","; echo "${gt[*]}");
#columns="${columns},${gt}";
}

#columns -> where filter is applied (keeping chrom,ref,alt to make sure loc same variants are pulled from dbs)
#use `gemini db_info <ensemble.db> | grep variants` to find out columns available in variants table
#there are 3 tables in gemini db: variants, variant_impacts, samples
columns="chrom,start,ref,alt,impact_severity,clinvar_pathogenic,clinvar_sig,gnomad_af_popmax,(gts).(*),(gt_depths).(*),(gt_alt_depths).(*)"
echo $columns | tr "," "\t"


#for each variant (location), 2 lines one from each db will be printed in consecutive lines
#if a variant was not found in a db, then a empty line will be returned for that line/db

while read -a l; do
gemini region --show-samples --reg "${l[0]}:${l[1]}-${l[2]}" --columns "$columns" $1 ;
gemini region --show-samples --reg "${l[0]}:${l[1]}-${l[2]}" --columns "$columns" $2 ;
echo;
done < $3