#!/bin/bash
#usage: compare_reports.sh <family.wes/clinical.csv> <family.wes/clinical.csv>




#https://www.biostars.org/p/347039/
db_extract () {
parallel  --colsep "\t"  'gemini region --reg "{1}:{2}-{3}" --header --columns "ref,alt,baseqranksum,dp,db,ac,impact_severity,clinvar_pathogenic,clinvar_sig"  $1' :::: $2 >> $3
#parallel --colsep "\t"  'gemini region --reg "{1}:{2}-{3}" --header --columns "ref,alt,baseqranksum,dp,db,ac,impact_severity,clinvar_pathogenic,clinvar_sig"  *-ensemble.db' :::: uniq.1.pos >> uniq.1.dblines
}

report_type=`basename $1 .csv | cut -d "." -f2`;
prefix1=`basename $1 .csv`;
prefix2=`basename $2 .csv`;

echo "given arguments ->";
echo "csv: $1 and $2";
echo "dbs: $3 and $3";
echo "report type = ${report_type}";

if [ "${report_type}" != "$(basename $2 .csv | cut -d "." -f2)" ]; then
	echo "report type (clinical or regular) should be same for both files. exiting!";
	exit;
fi;

grep -v "Position" $1 | cut -d "," -f1 | sort -k1,1 > 1.pos
grep -v "Position" $2 | cut -d "," -f1 | sort -k1,1 > 2.pos

left="${prefix1}.uniq.pos";
right="${prefix2}.uniq.pos";
common="${prefix1}.${prefix2}.common.pos";

comm -12 1.pos 2.pos | awk -vFS=":" -vOFS="\t" '{ split($1,a,"\""); print a[2],int($2)-1,int($2); }' > ${common}
comm -23 1.pos 2.pos | awk -vFS=":" -vOFS="\t" '{ split($1,a,"\""); print a[2],int($2)-1,int($2); }' > ${left}
comm -13 1.pos 2.pos | awk -vFS=":" -vOFS="\t" '{ split($1,a,"\""); print a[2],int($2)-1,int($2); }' > ${right}

echo "common pos : `wc -l ${common} | cut -d " " -f1`";
echo "pos only in $1: `wc -l ${left} | cut -d " " -f1`";
echo "pos only in $2: `wc -l ${right} |cut -d " " -f1`";