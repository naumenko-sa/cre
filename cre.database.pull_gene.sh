#!/bin/bash

gene=$1
database_prefix=$2
prefix=`date +%Y-%m-%d`

head -n1 $database_prefix.c4r.variant_wise.csv > $prefix.$gene.variant_wise.csv
cat $database_prefix.c4r.variant_wise.csv | grep $gene >> $prefix.$gene.variant_wise.csv

cat ~/cre/cre.database.header1 > $prefix.$gene.sample_wise.csv
cat $database_prefix.c4r.sample_wise.csv | grep $gene | awk -F "," '{print $NF","$0}' | sort >> $prefix.$gene.sample_wise.csv
