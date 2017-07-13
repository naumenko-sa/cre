#!/bin/bash

#package all reports to send as one archive

current_dir=`pwd`

destination=CHEO_`date +%F`

mkdir $destination

source_dir=/hpf/largeprojects/ccm_dccforge/dccforge/results

cd $source_dir

for dir in *x;do cd $dir; for f in *;do cp $f/$f.csv $current_dir/$destination;done;cd ..;done;

cd $current_dir

tar czf ${destination}.tar.gz -C $current_dir $destination

