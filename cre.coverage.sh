#!/bin/bash

#average coverage and % of bases covered at 10x
#total 35,405,680, 22387 genes

for f in *.coverage;do cat $f  | sed 1d | awk '{if ($9<10) print $0}'  | awk '{print $4}'  > $f.bad_genes.txt;done;
for f in *.bad_genes.txt; do for gene in `cat $f`;do cat protein_coding_genes.exons.fixed.bed | grep $gene >> $f.exon_lengths;done;done;
for f in *.exon_lengths;do echo $f;cat $f | awk '{sum+=$3-$2+1}END{print sum}';done;
