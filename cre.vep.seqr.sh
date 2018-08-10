#/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=5
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

#module load perl/5.20.1

reference=`which gatk-launch | sed s/"bin\/gatk-launch"/"\/genomes\/Hsapiens\/GRCh37"/`
echo ${reference}
bcbio_reference=`which gatk-launch | sed s/"bin\/gatk-launch"/"bcbio"/`
vep_reference=$(readlink -f `which vep`| sed s/"\/vep"//)

unset PERL5LIB && export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin:$PATH && /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin/vep --vcf -o stdout \
    -i $vcf --fork 16 --species homo_sapiens --no_stats --cache --offline --dir /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/vep --symbol --numbers --biotype --total_length \
    --canonical --gene_phenotype --ccds --uniprot --domains --regulatory --protein --tsl --appris --af --max_af --af_1kg --af_esp --af_gnomad --pubmed --variant_class \
    --allele_number \
    --fields Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,ALLELE_NUM,DISTANCE,STRAND,FLAGS,VARIANT_CLASS,SYMBOL_SOURCE,HGNC_ID,CANONICAL,TSL,APPRIS,CCDS,ENSP,SWISSPROT,TREMBL,UNIPARC,REFSEQ_MATCH,SOURCE,GIVEN_REF,USED_REF,BAM_EDIT,GENE_PHENO,SIFT,PolyPhen,DOMAINS,HGVS_OFFSET,AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,AA_AF,EA_AF,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,MAX_AF,MAX_AF_POPS,CLIN_SIG,SOMATIC,PHENO,PUBMED,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,LoF,LoF_filter,LoF_flags,LoF_info,Polyphen2_HVAR_pred,CADD_phred,SIFT_pred,FATHMM_pred,MutationTaster_pred,MetaSVM_pred \
    --fasta /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa.gz \
    --plugin LoF,human_ancestor_fa:/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/variation/human_ancestor.fa.gz,\
             loftee_path:/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/share/ensembl-vep-91.2-0 \
    --plugin dbNSFP,/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/variation/dbNSFP.txt.gz,Polyphen2_HVAR_pred,CADD_phred,SIFT_pred,FATHMM_pred,MutationTaster_pred,MetaSVM_pred \
    --plugin MaxEntScan,/hpf/largeprkjects/ccmbio/naumenko/tools/bcbio/anaconda/share/maxentscan-0_2004.04.21-0 \
    --plugin SpliceRegion --sift b --polyphen b --hgvs --shift_hgvs 1 --merged | \
     sed '/^#/! s/;;/;/g'  > `echo $vcf | sed s/vcf.gz/seqr.vcf.gz/`
