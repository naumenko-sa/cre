#!/bin/bash
# runs of HETEROZYGOUS VARIANTS MAF 5% variant in proband (CH0620) affected sib (CH0621) and in both
# filters DP=50, QUAL>=500
# gt_types 1 = HET
# usage:
# cre.roh.naive.sh sample gemini.db [maf=0.05]
# output:
# - sample.roh_variants.tsv - list of variants, print with tee for debugging
# - stdout: list of ROH

sample=$1
echo "chrom,pos,ref,alt,impact,qual,dp,gene,maf,gts."$sample",gt_types."$sample",gt_alt_depths."$sample",stretch_length_variants,stretch_length_bp,stretch_id,stretch_genes" | tee $sample.roh_variants.tsv

maf=0.05
if [ -n "$3" ]
then
    maf=$3
fi

gemini query -q "select chrom,start+1 as pos, ref, alt,impact,qual,depth, gene, max_aaf_all as maf, gts."$sample",gt_types."$sample",gt_alt_depths."$sample" from variants where
type='snp' and depth>=10 and qual>=100 and max_aaf_all<="$maf --gt-filter "gt_types."$sample" != 2" $2 | sed s/"\t"/","/g | sort -t "," -k1,1n -k2,2n \
| tee -a $sample.roh_variants.tsv | awk -F "," '
BEGIN{
    prev=1;
    prev_gene="";
    prev_chrom="";
    stretch_length_variants=0;
    stretch_length_bp=0;
    stretch_id="";
    stretch_genes="";
}
{
    genotype=$11;
    if(genotype != 1 || $1 != prev_chrom){
	stretch_length_variants=0;
	stretch_length_bp=0;
	stretch_id=0;
	stretch_genes="";
	prev_gene="";
    }else{
        if(prev==0){
	    stretch_length_variants=stretch_length_variants+1;
	    stretch_length_bp=$2-stretch_id+1;
	    if ($8 != prev_gene && $8 != ""){
    		stretch_genes=stretch_genes","$8;
    	    }
    	    prev_gene=$8;
	}else{
    	    stretch_length_variants=1;
    	    stretch_id=$2;
    	    stretch_length_bp=$2-stretch_id+1;
    	    stretch_genes=$8;
    	    prev_gene=$8;
	};
    }
    prev=genotype;
    prev_chrom=$1;
    print $0","stretch_length_variants","stretch_length_bp","stretch_id",\""stretch_genes"\"";
}' | grep -v "0$" | awk -F ',' '{ if ($13>=10) print $0;}' \
| sort -t "," -k1,1 -k15,15n \
| awk -F "," '
BEGIN{
    prev_chr="";
    prev_stretch_id="";
    prev_stretch="";
}
{
    if($1 != prev_chr){
	print prev_stretch;
    }else{
	if (prev_stretch_id != $15){
	    print prev_stretch;
	}
    }
    prev_stretch=$0;
    prev_chr=$1;
    prev_stretch_id=$15;
}' | grep -v "^$"

