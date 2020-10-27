#!/bin/bash

#usage: sh gemini_compare.sh <crg2-ensemble.db> <crg-ensemble.db> <pos.bed>

#$1 and $2 are ensemble dbs to extract variants
#$3 three-column(0-start) bed intervals (this could be common/uniq variants)
#$4-end any number of space-seperated column names to include in query 
	#(make sure column names are available in db; gemini outputs error for column not in db )


nargs=$#;
arr=(`echo $@ | awk '{ for(i=1;i<=NF;i++) print $i; }'`);

if [ $nargs -lt 3 ]; then
	echo "You must pass at least 3 arguments to the script."
	echo "usage:"
	echo "sh gemini_compare.sh ensemble1.db ensemble2.db pos.bed [space seperated column names in ensemble variants db]"
	echo "Exiting!"
	exit;
fi;

if [ $nargs -gt 3 ]; then 
	n=`expr $nargs - 1`;
	arg_columns=();
	for i in `seq 3 $n`; do
		arg_columns+=(${arr[$i]});
	done;
	arg_columns=","$(IFS=","; echo "${arg_columns[*]}");
else	
	arg_columns="";
fi;


#assuming the same sample naming in both the dbs use this instead of gts.(*)
get_samples () {
a=($(gemini query -q "select name from samples;" ${arr[0]}));
for i in ${a[*]}; do gt+=($(echo "gts.${i},gt_alt_depths.${i},gt_depths.${i}")); done;
gt=","$(IFS=","; echo "${gt[*]}");
}

#columns -> where filter is applied (keeping chrom,ref,alt to make sure loc same variants are pulled from dbs)
#use `gemini db_info <ensemble.db> | grep variants` to find out columns available in variants table
#there are 3 tables in gemini db: variants, variant_impacts, samples
#
fixed_columns="chrom,start,end,ref,alt";
filter_columns=",impact_severity,clinvar_pathogenic,clinvar_sig,gnomad_af_popmax";
get_samples #use ${gt} for genotype.sample columns
columns=${fixed_columns}${filter_columns}${arg_columns}${gt};
echo "${columns},db" | tr "," "\t"


#for each variant (location), 2 lines one from each db will be printed in consecutive lines
#if a variant was not found in a db, then a empty line will be returned for that line/db

while read -a l; do
ref=${l[3]};
alt=${l[4]};
echo -e "`gemini region  --reg "${l[0]}:${l[1]}-${l[2]}"  --columns "$columns"  --filter "ref='$ref' and alt='$alt'" ${arr[0]}` \t db1"
echo -e "`gemini region  --reg "${l[0]}:${l[1]}-${l[2]}"  --columns "$columns"  --filter "ref='$ref' and alt='$alt'" ${arr[1]}` \t db2"
echo;
done < ${arr[2]};
