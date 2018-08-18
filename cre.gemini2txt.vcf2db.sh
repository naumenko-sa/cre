#!/bin/bash
#  exports gemini.db database to gemini.db.txt file
#  database schema: https://gemini.readthedocs.io/en/latest/content/database_schema.html#the-variants-table
#  when using v.chr = g.chr AND v.gene = g.gene it becomes very slow
#  by default bcbio writes PASS only variants to the database

#  example call: cre.gemini2txt.sh S28-ensemble.db 5 ALL
#  when using vcfanno/vcfdb loader some fields are different
#  for some reason \n in the query string does not work here

#PBS -l walltime=1:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

if [ -z $file ]
then
    file=$1
fi

#10 reads for WES 5 reads for RNA-seq
depth_threshold=$2

severity_threshold=$3
#echo $severity_threshold

gemini query -q "select name from samples" $file > samples.txt

sQuery="select \
        variant_id as Variant_id,\
        ref as Ref,\
        alt as Alt,\
        impact as Variation,\
        dp as Depth,\
        qual as Quality,\
        gene as Gene,\
        ensembl_geneid as Ensembl_gene_id,\
        clinvar_sig as Clinvar,\
        transcript as Ensembl_transcript_id,\
        aa_length as AA_position,\
        exon as Exon,\
        domains as Protein_domains,\
        rs_ids as rsIDs,\
        af_1kg as Af_1000g,\
        gnomad_af_es as Gnomad_af_es,\
        gnomad_af_gs as Gnomad_af_gs,\
        max_af as Max_af,\
        gnomad_ac_gs as Gnomad_ac_gs,\
        gnomad_ac_es as Gnomad_ac_es,\
        gnomad_hom_gs as Gnomad_hom_gs,\
        gnomad_hom_es as Gnomad_hom_es,\
        sift_score as Sift_score,\
        polyphen_score as Polyphen_score,\
        cadd_phred as Cadd_score,\
        chrom as Chrom,\
        start+1 as Pos,\
        aa_change as AA_change,\
        hgvsc as Codon_change,\
        af_esp_aa as Evs_af_aa,\
        af_esp_ea as Evs_af_ea,\
        af_esp_all as Evs_af_all,\
        phylop20way_mammalian as Conserved_in_20_mammals,
        gts,"

while read sample;
do
	sQuery=$sQuery"gts."$sample","
	sQuery=$sQuery"gt_alt_depths."$sample","
	sQuery=$sQuery"gt_depths."$sample","
done < samples.txt

# gene_detailed may contain 2 records per single transcript - because of synonymous gene names, and some genes may have None in the name,for example TSRM
# https://groups.google.com/forum/#!topic/gemini-variation/U3uEvWCzuQo
# v.depth = 'None' see https://github.com/chapmanb/bcbio-nextgen/issues/1894

if [[ "$severity_threshold" == "ALL" ]]
then
#used for RNA-seq = 20k variants in the report
    severity_filter=""
#use for WES = 1k variants in the report
else
    severity_filter=" and impact_severity<>'LOW' "
fi

#max_aaf_all frequency is from gemini.conf and does not include gnomad WGS frequencing, gnomad WES only
#gnomad_af includes gnomad WGS
sQuery=$sQuery"hgvsc as Nucleotide_change_ensembl,\
		hgvsp as Protein_change_ensembl,\
        old_multiallelic as Old_multiallelic
		from variants \
        where \
	        (dp >= "$depth_threshold" or dp = '' or dp is null) "$severity_filter" and max_af < 0.01"

#echo $sQuery

gemini query --header -q "$sQuery" $file
