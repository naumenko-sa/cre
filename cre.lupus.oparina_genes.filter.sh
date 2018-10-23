#!/bin/bash

bname=`basename $1 .csv`

head -n1 $1 > $bname.lupus.csv
cat $1 | egrep "\"()\"" >> $bname.lupus.oparina.csv
#cat $1 | awk -F ',' '$6 ~/^(AP3B1|BLOC1S6|CD27|GATA2|ITK|LYST|NLRC4|PRF1|RAB27A|SH2D1A|SLC7A7|STX11|STXBP2|UNC13D|XIAP)$/{print}' > $bname.lupus.csv
