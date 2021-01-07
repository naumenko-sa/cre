#!/bin/bash

#usage: sh gemini_compare.sh <crg2-ensemble.db> <prefix1> <crg-ensemble.db> <prefix2> <pos.bed> [<space seperated column names>]

#$1 and $3 are ensemble dbs to extract variants
#$2 and $4 are prefixes appended as values to the retrieved db rows
#$5 three-column(0-start) bed intervals (this could be common/uniq variants)
#$6-end any number of space-seperated column names to include in query 
	#(make sure file names in the previous arguments don't have space as part of the names)
	#(tries to check if the column names are available in each db, if not uses "00" as column name) 

nargs=$#;
arr=(`echo $@ | awk '{ for(i=1;i<=NF;i++) print $i; }'`);

if [ $nargs -lt 5 ]; then
	echo "You must pass at least 3 arguments to the script."
	echo "usage:"
	echo "sh gemini_compare.sh ensemble1.db prefix 1 ensemble2.db prefix2 pos.bed [space seperated column names in ensemble variants db]"
	echo "Exiting!"
	exit;
fi;

if [ $nargs -gt 5 ]; then 
	n=`expr $nargs - 1`;
	arg_columns=();
	for i in `seq 5 $n`; do
		arg_columns+=(${arr[$i]});
	done;
	arg_str=","$(IFS=","; echo "${arg_columns[*]}");
else	
	arg_str="";
fi;

#assuming the same sample naming in both the dbs use this instead of gts.(*)
get_samples () {
a=($(gemini query -q "select name from samples;" ${arr[0]}));
for i in ${a[*]}; do gt+=($(echo "gts.${i},gt_alt_depths.${i},gt_depths.${i}")); done;
gt=","$(IFS=","; echo "${gt[*]}");
}

check_arg_columns () {

#pass db as first argument
for i in ${arg_columns[*]}; do 
	a=`gemini db_info $1 | grep -w "variants" | grep -w "$i"`;
	if [ ! -z "$a" ]; then  
		ret+=($i)
	else 
		ret+=("00")
	fi;
done;
if [ -z $ret ]; then 
	ret="";
else
	ret=","$(IFS=","; echo "${ret[*]}");
fi;
}

#columns -> where filter is applied (keeping chrom,ref,alt to make sure loc same variants are pulled from dbs)
#use `gemini db_info <ensemble.db> | grep variants` to find out columns available in variants table
#there are 3 tables in gemini db: variants, variant_impacts, samples
#

fixed_columns="chrom,start,end,ref,alt";
filter_columns=",impact_severity,clinvar_pathogenic,clinvar_sig,gnomad_af_popmax"; #gene,impact,ensembl_gene_id,callers";
column=${fixed_columns}${filter_columns};

get_samples #use ${gt} for genotype.sample columns

check_arg_columns ${arr[0]};
column1=${column}${ret}${gt}

unset ret;
check_arg_columns ${arr[2]};
column2=${column}${ret}${gt};

columns=${fixed_columns}${filter_columns}${arg_str}${gt}; 
# echo $column1;
# echo $column2;
# echo $columns;
echo "${columns},db" | tr "," "\t"


#for each variant (location), 2 lines one from each db will be printed in consecutive lines
#if a variant was not found in a db, then a empty line will be returned for that line/db

while read -a l; do
ref=${l[3]};
alt=${l[4]};
echo -e "`gemini region  --reg "${l[0]}:${l[1]}-${l[2]}"  --columns "${column1}"  --filter "ref='$ref' and alt='$alt'" ${arr[0]}` \t${arr[1]}"
echo -e "`gemini region  --reg "${l[0]}:${l[1]}-${l[2]}"  --columns "${column2}"  --filter "ref='$ref' and alt='$alt'" ${arr[2]}` \t${arr[3]}"
echo;
done < ${arr[4]};
