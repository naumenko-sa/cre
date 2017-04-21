# variant report generator

add_placeholder=function(variants,column_name,placeholder)
{
   variants[,column_name]=with(variants,placeholder)
   return(variants)
}

get_variants_from_file = function (filename)
{
    variants = read.delim(filename, stringsAsFactors=FALSE)
    return(variants)
}

# return Hom, Het, or - if == hom_reference
genotype2zygocity = function (genotype_str,ref)
{
      #genotype_str = "A|A|B"
      #genotype_str = "./." - call not possible
      #genotype_str = "TCA/."
      #genotype_str = "G"
      #genotype_str="A/A"
      #greedy
      genotype_str = gsub("|","/",genotype_str,fixed=T)
      genotype_str = gsub("./.","Insufficient_coverage",genotype_str,fixed=T)
      #genotype_str = gsub(".","NO_CALL",genotype_str,fixed=T)
      
      if(grepl("Insufficient_coverage",genotype_str)){
          result = genotype_str
      }else{
          ar = strsplit(genotype_str,"/",fixed=T)
          len = length(ar[[1]])
          if (len == 2)
          {
            if (ar[[1]][1] == ar[[1]][2]){
              if (ar[[1]][1] == ref)
                  result = "-"
              else
                  result = "Hom"
            }else
              result = "Het"
          }else{
            result = genotype_str
          }
      }
      return(result)
}

#output : family.ensemble.txt
create_report = function(family,samples)
{
    #test1: 3 samples in a family
    #family="166"
    #samples=c("166_3_5","166_4_10","166_4_8")
 
    #test2: 1 sample in a familty
    #family="NA12878-1"
    #samples=c("NA12878.1")
  
    #test3: 
    #setwd("~/Desktop/project_cheo/2017-01-30_dorin/1092R")
    #family="1092R"
    #samples = c("1092R_1613029","1092R_1613440","1092R_1613445")
    
    file=paste0(family,"-ensemble.db.txt")
  
    variants = get_variants_from_file(file)

    #Column1 - Position
    variants$Position=with(variants,paste(Chrom,Pos,sep=':'))
    
    #Column2 - UCSC link
    sUCSC1="=HYPERLINK(\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hgt.out3=10x&position="
    sUCSC2="\",\"UCSC_link\""
    variants$UCSC_Link=with(variants,paste(sUCSC1,Position,sUCSC2,")",sep=''))

    # Columns 3,4: Ref,Alt

    # Column 5 - Zygosity, column 7 - Burden
    # use new loader vcf2db.py - with flag  to load plain text
    # for genotype and depth - Noah
    # otherwise have to decode BLOB 
    # snappy decompression
    # https://github.com/arq5x/gemini/issues/700
    # https://github.com/lulyon/R-snappy
    for(sample in samples)
    {
        #DEBUG: gene = IL20RA
        #sample=samples[1]
        zygocity_column_name = paste0("Zygosity.",sample)
        #t = lapply(variants[,paste0("gts.",sample),"Ref"],genotype2zygocity)
        #t = lapply(variants[,paste0("gts.",sample),"Ref"],genotype2zygocity)
        t=unlist(mapply(genotype2zygocity,variants[,paste0("gts.",sample)],variants[,"Ref"]))
        variants[,zygocity_column_name] = unlist(t)
    
        burden_column_name = paste0("Burden.",sample)
        t = subset(variants, get(zygocity_column_name) == 'Hom' | get(zygocity_column_name) == 'Het',select=c("Ensembl_gene_id",zygocity_column_name))
        df_burden = count(t,'Ensembl_gene_id')    
        colnames(df_burden)[2] = burden_column_name
        variants = merge(variants,df_burden,all.x=T)
        variants[,burden_column_name][is.na(variants[,burden_column_name])] = 0
    }
    
    # Column6 - Gene
    
    # Column8 - gts
    
    # Column9 - Variation
    
    # Column 10 -  Info_ensembl
    variants = add_placeholder(variants,"Info_ensembl","Info_ensembl")
    impact_file=paste0(family,"-ensemble.db.impacts.txt")
    impacts = get_variants_from_file(impact_file)
    
    for (i in 1:nrow(variants))
    {
        #debug: i=1  
        v_id = variants[i,"Variant_id"]
        gene = variants[i,"Gene"]
        gene_impacts = subset(impacts,variant_id==v_id,select=c("exon","vep_hgvsc","vep_hgvsp"))
        gene_impacts$gene = rep(gene,nrow(gene_impacts))
        
        gene_impacts$exon[gene_impacts$exon=='']='NA'
        
        gene_impacts = gene_impacts[c("gene","exon","vep_hgvsc","vep_hgvsp")]
        
        v_impacts = paste0(gene_impacts$gene,":exon",gene_impacts$exon,":",gene_impacts$vep_hgvsc,":",gene_impacts$vep_hgvsp)
        s_impacts = paste(v_impacts,collapse=",")
      
        variants[i,"Info_ensembl"] = s_impacts
    }
    
    # Column 11 - Protein_change_ensembl
    
    # Column 12,13 - Info_refseq, Protein_change_refseq
    # refseq_impacts are pasted in merge_reports function
    variants = add_placeholder(variants,"Info_refseq","Info_refseq")
    variants = add_placeholder(variants,"Protein_change_refseq","NA")
    
    # Columns 14,15 - Depth, Qual_depth

    # Column 16 - Alt_depth - from v.gt_alt_depth
    # when multiple callers used, AD is not set and fixed in merge_reports function
    for(sample in samples)
    {
        new_name = paste0("Alt_depths.",sample)
        setnames(variants, paste0("gt_alt_depths.",sample),new_name)
    }

    # Column 17 - Trio_coverage - fixed in merge_reports function
    variants = add_placeholder(variants,"Trio_coverage","")
    n_sample = 1
    prefix = ""
    
    #order gts column in the same way as in samples
    variants$gts=""
    for(sample in samples)
    {
        column = paste0("gt_depths.",sample)
        
        if (n_sample>1) prefix="/"
        
        variants$Trio_coverage = with(variants,paste0(Trio_coverage,prefix,get(column)))
      
        column = paste0("gts.",sample)
        
        if (n_sample>1) prefix=","
        
        variants$gts = with(variants,paste0(gts,prefix,get(column)))
      
        n_sample = n_sample+1
    }
    
    # Column 18 - Ensembl_gene_id

    # Column 19 - Gene_description
    gene_descriptions = read.delim2(paste0(default_tables_path,"/ensembl_w_description.txt"), stringsAsFactors=F)
    variants = merge(variants,gene_descriptions,by.x = "Ensembl_gene_id",by.y = "ensembl_gene_id",all.x=T)
    
    # Column 20 - Omim_gene_description - from omim text file
    # omim.forannotation2 previously
    # I use my copy of omim first, then look into cre, because I don't want to distribute omim
    # through my github.
    omim_file_name = paste0(default_tables_path,"/omim.txt")
    omim_file_name_local = paste0(reference_tables_path,"/omim.txt")
    
    if (file.exists(omim_file_name_local)) omim_file_name = omim_file_name_local
  
    omim = read.delim2(omim_file_name_local, stringsAsFactors=F)
    variants = merge(variants,omim,all.x=T)

    # Column 21 - Omim_inheritance 
    omim_file_name = paste0(default_tables_path,"/omim_inheritance.txt")
    omim_file_name_local = paste0(reference_tables_path,"/omim_inheritance.txt")
    
    if (file.exists(omim_file_name_local)) omim_file_name = omim_file_name_local
    
    omim_inheritance = read.csv(omim_file_name, sep="",stringsAsFactors = F)
    variants = merge(variants,omim_inheritance,all.x=T)

    # Column 22 - Orphanet
    # previous name - orphanet.deduplicated.txt
    orphanet_file_name = paste0(default_tables_path,"/orphanet.txt")
    orphanet_file_name_local = paste0(reference_tables_path,"/orphanet.txt")
    
    if (file.exists(orphanet_file_name_local)) orphanet_file_name = orphanet_file_name_local
    
    orphanet = read.delim(orphanet_file_name, stringsAsFactors=F)  
    variants = merge(variants,orphanet,all.x=T)
    
    # Column 23 - Clinvar
    
    # Column 24 - Ensembl_transcript_id
    
    # Column 25 - AA_position
    
    # Column 26 - Exon
    
    # Column 27 - Pfam_domain
    
    # Column 28, 29 = Frequency_in_C4R, Seen_in_C4R_samples
    variants = add_placeholder(variants,"Frequency_in_C4R","Frequency_in_C4R")
    variants = add_placeholder(variants,"Seen_in_C4R_samples","Seen_in_C4R_samples")

    # Column 30 - rsIds

    # Columns 31-36 - population frequencies

    # Columns 37-38, Exac scores
    exac_scores_file = paste0(default_tables_path,"/exac_scores.txt")
    exac_scores = read.delim(exac_scores_file, stringsAsFactors=F)
    variants = merge(variants,exac_scores,all.x=T)

    # Column 39 - Exac_het
    # Column 40 - Exac_hom_alt
    
    # Column 41 - Conserved in 29 mammals instead of phastcons
    #https://www.biostars.org/p/150152/

    # Column 42-43-44: sift,polyphen,cadd scores

    
    # Columns 45,46 - imprinting
    imprinting_file_name = paste0(default_tables_path,"/imprinting.txt")
    imprinting = read.delim(imprinting_file_name, stringsAsFactors=F)
    variants = merge(variants,imprinting,all.x=T)
    
    # Column 47 - pseudoautosomal
    pseudoautosomal_file_name = paste0(default_tables_path,"/pseudoautosomal.txt")
    pseudoautosomal = read.delim(pseudoautosomal_file_name, stringsAsFactors=F)
    variants = merge(variants,pseudoautosomal,all.x=T)
    
    # replace -1 with 0
    for (field in c("EVS_maf_aa","EVS_maf_ea","EVS_maf_all","Maf_1000g","Exac_maf","Maf_all","Exac_het","Exac_hom_alt","Trio_coverage"))
    {
        variants[,field] = with(variants,gsub("-1","0",get(field),fixed=T))  
    }

    for (field in c(paste0("Alt_depths.",samples)))
    {
        variants[,field] = with(variants,gsub("-1",NA,get(field),fixed=T))  
    }

    select_and_write2(variants,samples,paste0(family,".create_report"))
}

# column selection and order
select_and_write = function(variants,samples,prefix)
{
    variants = variants[c(c("Position","UCSC_Link","Ref","Alt"),
                        paste0("Zygosity.",samples),c("Gene"),
                        paste0("Burden.",samples),c("gts","Variation","Info_ensembl","Protein_change_ensembl","Info_refseq","Protein_change_refseq","Depth","Quality"),
                        paste0("Alt_depths.",samples),
                        c("Trio_coverage","Ensembl_gene_id","Gene_description","Omim_gene_description","Omim_inheritance",
                          "Orphanet", "Clinvar","Ensembl_transcript_id","AA_position","Exon","Pfam_domain",
                          "Frequency_in_C4R","Seen_in_C4R_samples","rsIDs","Maf_1000g","EVS_maf_aa","EVS_maf_ea","EVS_maf_all",
                          "Exac_maf","Maf_all", "Exac_pLi_score","Exac_missense_score","Exac_het","Exac_hom_alt",
                          "Conserved_in_29_mammals","Sift_score","Polyphen_score","Cadd_score",
                          "Imprinting_status","Imprinting_expressed_allele","Pseudoautosomal"))]
  
    write.table(variants,paste0(prefix,".txt"),quote=F,sep = ";",row.names=F)  
}

# writes in CSV format
select_and_write2 = function(variants,samples,prefix)
{
  variants = variants[c(c("Position","UCSC_Link","Ref","Alt"),paste0("Zygosity.",samples),c("Gene"),
                        paste0("Burden.",samples),c("gts","Variation","Info_ensembl","Protein_change_ensembl","Info_refseq","Protein_change_refseq","Depth","Quality"),
                        paste0("Alt_depths.",samples),c("Trio_coverage","Ensembl_gene_id","Gene_description","Omim_gene_description","Omim_inheritance",
                                                        "Orphanet", "Clinvar","Ensembl_transcript_id","AA_position","Exon","Pfam_domain",
                                                        "Frequency_in_C4R","Seen_in_C4R_samples","rsIDs","Maf_1000g","EVS_maf_aa","EVS_maf_ea","EVS_maf_all",
                                                        "Exac_maf","Maf_all", "Exac_pLi_score","Exac_missense_score","Exac_het","Exac_hom_alt",
                                                        "Conserved_in_29_mammals","Sift_score","Polyphen_score","Cadd_score",
                                                        "Imprinting_status","Imprinting_expressed_allele","Pseudoautosomal"))]
  
  write.csv(variants,paste0(prefix,".csv"),row.names = F)  
}

fix_column_name = function(column_name)
{
    if(grepl("^[0-9]",column_name))
    {
      column_name = paste0("X",column_name)
    }
    return(column_name)
}

# merges ensembl, gatk-haplotype reports
merge_reports = function(family,samples)
{
    # test:
    # setwd("/home/sergey/Desktop/project_cheo/2016-11-09_rerun10")
    # family = "166"
    # mind the samples order: it will influence the Trio
    # samples=c("166_3_5","166_4_10","166_4_8")
  
    ensemble_file = paste0(family,".create_report.csv")
    
    ensemble = read.csv(ensemble_file, stringsAsFactors=F)
    ensemble$superindex=with(ensemble,paste(Position,Ref,Alt,sep='-'))
    
    refseq_file = paste0(family,".refseq.txt")
    refseq = read.delim(refseq_file, stringsAsFactors=F,na.strings = "")
    ensemble = merge(ensemble,refseq,by.x = "superindex", by.y="superindex",all.x = T)
    
    for (i in 1:nrow(ensemble))
    {
        if (is.na(ensemble[i,"Info_refseq_no_gene"]))
        {
              ensemble[i,"Info_refseq"] = NA
              ensemble[i,"Protein_change_refseq"] = NA
        }
        else
        {
            v_impacts = strsplit(ensemble[i,"Info_refseq_no_gene"],",",fixed=T)[[1]]
            gene = ensemble[i,"Gene"]
            ensemble[i,"Info_refseq"]=paste(paste(gene,v_impacts,sep=":"),collapse=",")
            for (impact in v_impacts)
            {
                if (grepl(":NP_",impact,fixed = T))
                {
                    v_subimpacts = strsplit(impact,":",fixed=T)[[1]]
                    ensemble[i,"Protein_change_refseq"] = v_subimpacts[5]
                    break
                }
            }
        }
    }
    
    gatk_file = paste0(family,"-gatk-haplotype-annotated-decomposed.table")
    gatk = read.delim(gatk_file, stringsAsFactors=F)
    gatk$superindex=with(gatk,paste(paste0("chr",CHROM,":",POS),REF,ALT,sep='-'))
    gatk[c("CHROM","POS","REF","ALT")]=NULL
    
    ensemble = merge(ensemble,gatk,by.x = "superindex", by.y="superindex",all.x = T)
    
    ensemble$Depth = ensemble$DP
    n_sample = 1
    prefix = ""
    ensemble$Trio_coverage=""
    
    for(sample in samples)
    {
        #R fixes numerical column names with X?
        #what if sample is not numerical
      
        column = fix_column_name(sample)  
        column = paste0(column,".DP")
        
        if (n_sample>1) prefix="/"
        
        ensemble$Trio_coverage = with(ensemble,paste0(Trio_coverage,prefix,get(column)))
      
        column = paste0("Alt_depths.",sample)
        column_gatk = fix_column_name(sample)
        column_gatk = paste0(column_gatk,".AD")
        
        ensemble[,column] = ensemble[,column_gatk]
      
        n_sample = n_sample+1
    }
    
    for (i in 1:nrow(ensemble))
    {
        for (sample in samples)
        {
            field = paste0("Alt_depths.",sample)
            
            ensemble[i,field]=strsplit(ensemble[i,field],",",fixed=T)[[1]][2]
        }
    }
    
    for (sample in samples)
    {
        ensemble[c("DP",paste0(fix_column_name(sample),".DP"),paste0(fix_column_name(sample),".AD"))]=NULL
    }
    
    freebayes_file = paste0(family,"-freebayes-annotated-decomposed.table")
    freebayes = read.delim(freebayes_file, stringsAsFactors=F)
    freebayes$superindex=with(freebayes,paste(paste0("chr",CHROM,":",POS),REF,ALT,sep='-'))
    freebayes[c("CHROM","POS","REF","ALT")]=NULL
    ensemble = merge(ensemble,freebayes,by.x = "superindex", by.y="superindex",all.x = T)
    
    for (i in 1:nrow(ensemble))
    {
        
        #if(grepl("NA",ensemble[i,"Trio_coverage"]))
        #wrong: a variant may be called by gatk with 10/10/NA,
        #and freebayes will destroy coverage info
        if (str_count(ensemble[i,"Trio_coverage"],"NA") == length(samples))
        {
            ensemble[i,"Depth"] = ensemble[i,"DP"]
            for (sample in samples)
            {
                field_depth = paste0("Alt_depths.",sample)
                field_bayes = paste0(fix_column_name(sample),".AO")
                #field_bayes = paste0(sample,".AO")
                      
                ensemble[i,field_depth] = ensemble[i,field_bayes]
        
            }
            n_sample = 1
            prefix = ""
            ensemble[i,"Trio_coverage"]=""
            
            for(sample in samples)
            {
                column = paste0(fix_column_name(sample),".DP")
                if (n_sample>1) prefix="/"
                ensemble[i,"Trio_coverage"] = paste(ensemble[i,"Trio_coverage"],ensemble[i,column],sep = prefix)
              
                n_sample = n_sample+1
            }
        }
    }
    
    for (sample in samples)
    {
        ensemble[c("DP",paste0(fix_column_name(sample),".DP"),paste0(fix_column_name(sample),".AO"))]=NULL
    }
    
    platypus_file = paste0(family,"-platypus-annotated-decomposed.table")
    platypus = read.delim(platypus_file, stringsAsFactors=F)
    platypus$superindex=with(platypus,paste(paste0("chr",CHROM,":",POS),REF,ALT,sep='-'))
    platypus[c("CHROM","POS","REF","ALT")]=NULL
    ensemble = merge(ensemble,platypus,by.x = "superindex", by.y="superindex",all.x = T)
    
    for (i in 1:nrow(ensemble))
    {
      #if(grepl("NA",ensemble[i,"Trio_coverage"])) - wrong, may be 10/10/NA in gatk
      #if (ensemble[i,"Trio_coverage"]=="NA/NA/NA")
      if (str_count(ensemble[i,"Trio_coverage"],"NA") == length(samples))
      {
        ensemble[i,"Depth"] = ensemble[i,"TC"]
        for (sample in samples)
        {
          field_depth = paste0("Alt_depths.",sample)
          field_bayes = paste0(fix_column_name(sample),".NV")
          
          #sometimes freebayes has 10,10,10 for decomposed alleles
          ensemble[i,field_depth] = strsplit(ensemble[i,field_bayes],",",fixed=T)[[1]][1]
          
        }
        n_sample = 1
        prefix = ""
        ensemble[i,"Trio_coverage"]=""
        
        for(sample in samples)
        {
          column = paste0(fix_column_name(sample),".NR")
          if (n_sample>1) prefix="/"
          #sometimes freebayes has 10,10,10 for decomposed alleles
          cov_value=strsplit(ensemble[i,column],",",fixed=T)[[1]][1]
          ensemble[i,"Trio_coverage"] = paste(ensemble[i,"Trio_coverage"],cov_value,sep = prefix)
          
          n_sample = n_sample+1
        }
      }
    }
    
    for (sample in samples)
    {
        ensemble[c("TC",paste0(fix_column_name(sample),".NV"),paste0(fix_column_name(sample),".NR"))]=NULL
    }
    
    ensemble[,"Trio_coverage"] = with(ensemble,gsub("NA","0",get("Trio_coverage"),fixed=T))  
   
    for (i in 1:nrow(ensemble))
    {
        if (is.na(ensemble[i,"Depth"]))
        {
            l=strsplit(ensemble[i,"Trio_coverage"],"/")[[1]]
            ensemble[i,"Depth"]=sum(as.integer(l))
        }
        for (sample in samples)
        {
            field_depth = paste0("Alt_depths.",sample)
            if (is.na(ensemble[i,field_depth]))
                ensemble[i,field_depth]=0
        }
    }
    
    select_and_write2(ensemble,samples,paste0(family,".merge_reports"))
}

annotate_w_care4rare = function(family,samples)
{
    variants = read.csv(paste0(family,".merge_reports.csv"), stringsAsFactors=F)
  
    variants$superindex=with(variants,paste(Position,Ref,Alt,sep='-'))
    
    if(exists("seen_in_c4r_counts"))
    {
        variants = merge(variants,seen_in_c4r_counts,by.x = "superindex", by.y="superindex",all.x = T)
        variants$Frequency_in_C4R = variants$counts
        variants$counts=NULL
    }
    
    if(exists("seen_in_c4r_samples"))
    {
        variants = merge(variants,seen_in_c4r_samples,by.x = "superindex", by.y="superindex",all.x = T)
        variants$Seen_in_C4R_samples=variants$samples
    }
    
    select_and_write2(variants,samples,family)
}

library(stringr)
library(data.table)
library(plyr)

default_tables_path="~/cre"
reference_tables_path = "~/Desktop/reference_tables"

#load c4r information
seen_in_c4r_counts.txt = paste0(reference_tables_path,"/seen_in_c4r_counts.txt")
if (file.exists(seen_in_c4r_counts.txt))
{
    seen_in_c4r_counts = read.delim(seen_in_c4r_counts.txt, stringsAsFactors=F)
}
seen_in_c4r_samples.txt = paste0(reference_tables_path,"/seen_in_c4r_samples.txt")
if (file.exists(seen_in_c4r_samples.txt))
{
    seen_in_c4r_samples = read.csv(seen_in_c4r_samples.txt, stringsAsFactors=F, sep=";")
}

# R substitutes "-" with "." in sample names in columns so fix this in samples.txt
# sample names starting with letters should be prefixed by X in *.table
# for correct processing. most of them start with numbers, and R adds X automatically

args = commandArgs(trailingOnly = T)
family = args[1]

# DEBUG - replace with Ashkenazim trio
# setwd("~/Desktop/project_cheo/2017-03-16_NextSeq/")
# family="CHEO_0001"

setwd(family)

samples = unlist(read.table("samples.txt", stringsAsFactors=F))
samples = gsub("-",".",samples)
    
create_report(family,samples)
merge_reports(family,samples)
annotate_w_care4rare(family,samples)

setwd("..")
