args = commandArgs(trailingOnly = T)

#input_report = "182208.wes.2018-11-30.csv"

input_report = args[1]

variants = read.csv(input_report,stringsAsFactors = F)

lupus = read.csv("~/cre/data/lupus.csv",stringsAsFactors = F)

lupus.langefeld.gwas = read.csv("~/cre/data/lupus.langefeld.gwas.csv",stringsAsFactors = F)

lupus.eastasianssle.gwas = read.csv("~/cre/data/lupus.eastasianssle.gwas.csv", stringsAsFactors = F)

primary_immunodeficiency = read.csv("~/cre/data/primary_immunodeficiency.csv", stringsAsFactors = F)

periodic_fever = read.csv("~/cre/data/periodic_fever_syndromes.csv", stringsAsFactors = F)

variants$Lupus_panel = ifelse(variants$Ensembl_gene_id %in% lupus$ensembl_gene_id,"Lupus_panel",NA)
variants$Lupus_langefeld_gwas = ifelse(variants$Ensembl_gene_id %in% lupus.langefeld.gwas$ensembl_gene_id,"Lupus_langefeld_gwas",NA)
variants$Lupus_eastasianssle_gwas = ifelse(variants$Ensembl_gene_id %in% lupus.eastasianssle.gwas$ensembl_gene_id,"Lupus_eastasianssle_gwas",NA)
variants$Primary_immunodeficiency_panel = ifelse(variants$Ensembl_gene_id %in% primary_immunodeficiency$PanelAPP.EnsemblId.GRch37,"Primary_immunodeficiency_panel",NA)
variants$Periodic_fever_panel = ifelse(variants$Ensembl_gene_id %in% periodic_fever$PanelAPP.EnsemblId.GRch37,"Periodic_fever",NA)

output_report = sub(".csv",".lupus.csv",input_report)

write.csv(variants,output_report,row.names = F)