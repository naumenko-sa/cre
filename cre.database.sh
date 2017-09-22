#!/bin/bash

#$1 - input dir
#$2 - output dir

cd $1

for directory in *x
do 
    echo $directory
    cd $directory
    
    for family in *
    do
	cd $family
	cre.database.py $family
	mv *.c4r $2
	cd ..
    done
    cd ..
done

cd $2

prefix=`date +%Y-%m-%d`
cat ~/cre/cre.database.header *.c4r > $prefix.c4r.sample_wise.csv
cre.database_merge.py $prefix.c4r.sample_wise.csv $prefix.c4r.variant_wise.csv
#rm *.c4r

#create files for report generation
cat $prefix.c4r.variant_wise.csv | awk -F ',' '{print $1"-"$2"-"$3"\t"$(NF-1)}'  > ~/Desktop/reference_tables/seen_in_c4r_counts.txt
cat $prefix.c4r.variant_wise.csv | awk -F ',' '{print $1"-"$2"-"$3"\t"$NF}'  > ~/Desktop/reference_tables/seen_in_c4r_samples.txt

