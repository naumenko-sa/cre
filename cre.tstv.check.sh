#!/bin/bash

for f in *x/*;
do 
    if [ -f  ${f}/multiqc/multiqc_data/multiqc_data_final.json ]; 
    then 
	for g in `cat $f/multiqc/multiqc_data/multiqc_data_final.json | grep generalstats-tstv | awk -F ':' '{print $2}' | grep -v mqc | sed s/" "// | sed s/","//`;
	do 
	    echo $f $g >> tstv_check.txt;
	done;
    else
	for g in `cat $f/multiqc/multiqc_data/multiqc_bcftools_stats.txt | sed 1d | awk '{print $6}'`;
	do
	    echo $f $g  >> tstv_check.txt;
	done;
    fi
done;

cat tstv_check.txt | sort -n -k2,2 > tstv_check_sorted.txt