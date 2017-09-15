#!/bin/bash

gene=$1
database_prefix=$2
prefix=`date +%Y-%m-%d`

head -n1 $database_prefix.c4r.variant_wise.csv > $prefix.$gene.variant_wise
cat $database_prefix.c4r.variant_wise.csv | grep $gene >> $prefix.$gene.variant_wise

cat ~/cre/cre.database.header1 > $prefix.$gene.sample_wise
cat $database_prefix.c4r.sample_wise.csv | grep $gene | awk -F "," '{print $12","$0}' | sort >> $prefix.$gene.sample_wise
