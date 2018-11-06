#!/bin/bash

# prepares orpanet database for cre report generator
# result is orphanet.txt - should be 3602 genes

read_dom () {
    local IFS=\>
    read -d \< ENTITY CONTENT
}

wget -c http://www.orphadata.org/data/xml/en_product6.xml

while read_dom; 
do
    echo "$ENTITY => $CONTENT"
done < en_product6.xml > en_product6.xml.txt

cat en_product6.xml.txt | egrep -a "(Disorder id)|(Name lang)|ENSG" > en_product6.xml.parsed

cat en_product6.xml.parsed | awk '{if($0~"Disorder") {print $0;getline;print $0;}; if ($0~"ENSG") print $0;}' | grep -a -v "Disorder id" | awk -F '=> ' '{print $2}' > en_product6.xml.final

cat en_product6.xml.final | awk '{if($0 ~ "ENSG") {print $0"\t"dis}else{dis=$0}}' > orphanet.tmp

#there are cases of many diseases for single gene
cat orphanet.tmp | sort -k1,1 > orphanet.sorted.txt
cat orphanet.sorted.txt  | awk -F "\t" 'BEGIN{prev_gene="Ensembl_gene_id\tOrphanet";buf=""}{if(prev_gene != $1){print prev_gene"\t"buf;buf=$2;prev_gene=$1}else{buf=buf","$2;}}END{print prev_gene","buf'} > orphanet.txt

rm en_product6.* orphanet.tmp orphanet.sorted.txt
