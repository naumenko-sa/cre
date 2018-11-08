# get gnomad gene constraint scores
# - Gnomad_oe_lof_score
# - Gnomad_oe_mis_score
# output: gnomad_scores.csv
# run: Rscript ~/cre/cre.gnomad_scores.R
# https://macarthurlab.org/2018/10/17/gnomad-v2-1/

# install.packages("R.utils")
# bash:
# cd
# git clone https://github.com/naumenko-sa/bioscripts

source("~/bioscripts/genes.R")
library("R.utils")

gnomad_scores_url = "https://storage.googleapis.com/gnomad-public/release/2.1/ht/constraint/constraint.txt.bgz"
download.file(gnomad_scores_url,"gnomad_scores.txt.bgz")
gunzip("gnomad_scores.txt.bgz","gnomad_scores.txt")
gnomad_scores = read.delim("gnomad_scores.txt", stringsAsFactors=F)
gnomad_scores = gnomad_scores[,c("gene","transcript","canonical","oe_lof","oe_mis")]
gnomad_scores = gnomad_scores[gnomad_scores$canonical == "true",]

#still has a few duplicates
#gnomad_scores[duplicated(gnomad_scores$gene),]
gnomad_scores = gnomad_scores[!duplicated(gnomad_scores$gene),]

mart = init_mart_human()
get_protein_coding_genes(mart)

genes_transcripts = read.csv("genes.transcripts.csv", stringsAsFactors = F)
gnomad_scores = merge(gnomad_scores,genes_transcripts,by.x="transcript",by.y="Ensembl_transcript_id",all.x=T,all.y=F)

#some genes are absent in grch37
#gnomad_scores[is.na(gnomad_scores$Ensembl_gene_id),]

gnomad_scores = gnomad_scores[!is.na(gnomad_scores$Ensembl_gene_id),]
gnomad_scores = gnomad_scores[,c("Ensembl_gene_id","oe_lof","oe_mis")]

colnames(gnomad_scores) = c("Ensembl_gene_id","Gnomad_oe_lof_score","Gnomad_oe_mis_score")
write.csv(gnomad_scores,"gnomad_scores.csv",row.names = F)

file.remove("genes.transcripts.csv")
file.remove("gnomad_scores.txt")
