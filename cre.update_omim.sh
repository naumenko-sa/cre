# Updates OMIM annotations for the cre pipeline
# Use: navigate to omim directory and run ./update_omim.sh
# Input:
#  * mimTitles.txt genemap2.txt mim2gene.txt
# Output: 
#  * hgnc_$DATE.tsv
#    * a table used to map gene names with ensembl gene id's
#    * required file for cre.sh
#  * `prefix_hgnc_join_omim_phenos_yyyy-mm-dd.tsv`
#    * a table mapping each OMIM gene to all of its known HGNC gene symbol/previous symbols/synonyms
#    * allows for simple left-join to 'variant' column of CRE reports
#  * `prefix_remaining_unmapped_omim_with_info_yyyy-mm-dd.tsv`
#    * a list of the the unmapped OMIM genes

DATE=`date +%Y-%m-%d`
HGNC="hgnc_${DATE}.tsv"
OMIM_PREFIX="OMIM"

# get latest mappings from hgnc
# http params include default checkboxes, with ensembl gene id from hgnc, and ncbi gene id
wget \
"https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_pub_eg_id&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit" \
-O ${HGNC}

Rscript ~/cre/map_omim_hgnc.R ${HGNC} mimTitles.txt genemap2.txt mim2gene.txt ${OMIM_PREFIX}
