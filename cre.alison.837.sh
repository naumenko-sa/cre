#!/bin/bash

# get stretches of identical genotypes in trio for Alison

echo chrom,pos,ref,alt,gene,837_ED0054,837_ED0055,837_ED0056,polymorphism,stretch > 837.chr12.stretch.csv

gemini query \
-q "select chrom,start+1 as pos,ref,alt,gts.837_ED0054,gts.837_ED0055, gts.837_ED0056,gene from variants where chrom=12 and dp>20 and type='snp'" 837-ensemble.db \
| sed s/"G\/A"/"A\/G"/g | sed s/"T\/C"/"C\/T"/g | sed s/"G\/C"/"C\/G"/g | sed s/"C\/A"/"A\/C"/g | sed s/"T\/G"/"G\/T"/g | sed s/"T\/A"/"A\/T"/g \
| grep -v "\." | sort -n -k1,1 -k2,2 \
| awk '{if(($6==$7) && ($7==$8)){print $0"\tI"}else{print $0"\tP"}}' \
| awk 'BEGIN{prev="P";stretch=0}{if($9==prev && $9=="I"){stretch=stretch+1}else{stretch=0};prev=$9;print $0"\t"stretch}' | sed s/"\t"/","/g  >> 837.chr12.stretch.csv