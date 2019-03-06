#!/bin/bash
# runs of HOM MAF 5% variant in proband (CH0620) affected sib (CH0621) and in both
# gt_types=0 or 3 = HOM, 1 = HET
# usage:
# cre.roh.naive.sh sample gemini.db [maf=0.05]

sample=$1
echo "chrom,pos,ref,alt,impact,qual,dp,gene,maf,gts."$sample",gt_types."$sample",stretch_length_variants,stretch_length_bp,stretch_id,stretch_genes"

maf=0.05
if [ -n $3 ]
then
    maf=$3
fi

gemini query -q "select chrom,start+1 as pos, ref, alt,impact,qual,dp, gene, gnomad_af_popmax as maf, gts."$sample",gt_types."$sample" from variants where
type='snp' and dp>20 and qual>100 and gnomad_af_popmax<"$maf --gt-filter "gt_types."$sample" != 2" $2 | sort -k1,1n -k2,2n \
| awk -F "\t" '
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
    if($11 == 3){
	genotype=0
    }else{
	genotype=$11
    };
    if(genotype==1 || $1 != prev_chrom){
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
    print $0"\t"stretch_length_variants"\t"stretch_length_bp"\t"stretch_id"\t\""stretch_genes"\"";
}' | sed s/"\t"/","/g | grep -v "0$" | awk -F ',' '{ if ($12>9) print $0;}' \
| sort -k1,1 -k14,14n \
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
	if (prev_stretch_id != $14){
	    print prev_stretch;
	}
    }
    prev_stretch=$0;
    prev_chr=$1;
    prev_stretch_id=$14;
}' | grep -v "^$"

