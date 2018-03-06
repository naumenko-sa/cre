#!/bin/bash

#for s in `cat samples_ready.txt`;do if grep -q $s samples_w_vcf_consent.txt;then echo yes; else echo no;fi;done;
#for s in `cat samples_ready.txt`;do if grep -q $s samples_w_vcf_consent.txt;then echo $s >> 2upload.txt;fi;done;
#for s in `cat 2upload.txt`;do fam=`echo $s | awk -F '_' '{print $1}'`;pc=`cat samples_w_vcf_consent.txt | grep $s | awk '{print $1}'`;if [ -d $fam ]; then cd $fam;if [ -f $s.vcf ];then pc.upload_vcf.sh $pc $s.vcf;fi;cd ..;fi;done;

module load python/3.5.2
python3 ~/cre/pc.upload_vcf.py $1 $2  ~/work/project_cheo/pc.txt phenomecentral.org
