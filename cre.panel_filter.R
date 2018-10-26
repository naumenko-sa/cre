setwd("~/Desktop/work")

variants = read.csv("ivakine01.wgs.coding.2018-10-26.csv",stringsAsFactors = F)

panel = read.delim("Periodic fever syndromes.tsv", stringsAsFactors=F)

panel = panel[,c("Gene.Symbol","Model_Of_Inheritance","Phenotypes","UserRatings_Green_amber_red","EnsemblId.GRch37.")]

colnames(panel)=c("PanelAPP.Gene.Symbol","PanelAPP.Model_Of_Inheritance","PanelAPP.Phenotypes",
                  "PanelAPP.UserRatings_Green_amber_red", "PanelAPP.EnsemblId.GRch37")


variants.ens = variants[variants$Ensembl_gene_id %in% panel$PanelAPP.EnsemblId.GRch37,]
variants.gene = variants[variants$Gene %in% panel$PanelAPP.Gene.Symbol,]
