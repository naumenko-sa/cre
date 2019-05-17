# variant report generator
# Rscript ~/cre/cre.vcf2.db.R <family> noncoding|default=NULL,coding
# description of columns: 
# https://docs.google.com/document/d/1zL4QoINtkUd15a0AK4WzxXoTWp2MRcuQ9l_P9-xSlS4/edit?usp=sharing
add_placeholder <- function(variants, column_name, placeholder){
    variants[,column_name] <- with(variants, placeholder)
    return(variants)
}

# returns Hom / Het / - (for HOM reference)
genotype2zygosity <- function (genotype_str, ref){
    # test
    # genotype_str = "A|A|B"
    # genotype_str = "./." - call not possible
    # genotype_str = "TCA/."
    # genotype_str = "G"
    # genotype_str = "A/A"
    # greedy
    genotype_str <- gsub("|", "/", genotype_str, fixed = T)
    genotype_str <- gsub("./.", "Insufficient_coverage", genotype_str, fixed = T)
    
    if(grepl("Insufficient_coverage", genotype_str)){
      result <- genotype_str
    }else{
        ar <- strsplit(genotype_str, "/", fixed = T)
        len <- length(ar[[1]])
        if (len == 2){
            if (ar[[1]][1] == ar[[1]][2]){
                if (ar[[1]][1] == ref)
                    result <- "-"
                else
                    result <- "Hom"
            }else result <- "Het"
        }else result <- genotype_str
    }
    return(result)
}

# output : family.ensemble.txt
create_report <- function(family, samples){
    file <- paste0(family, ".variants.txt")
    variants <- read_delim(file, delim = "\t")
    
    impact_file <- paste0(family, ".variant_impacts.txt")
    impacts <- read_delim(impact_file, delim = "\t")
    
    #Column1 - Position
    variants <- variants %>% mutate(Position = paste(Chrom, Pos, sep = ':'))
    
    #Column2 - UCSC link
    sUCSC1 <- "=HYPERLINK(\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hgt.out3=10x&position="
    sUCSC2 <- "\",\"UCSC_link\")"
    
    variants <- variants %>% mutate(UCSC_Link = paste(sUCSC1, Position, sUCSC2, sep = ''))

    # Column3 = GNOMAD_Link
    sGNOMAD1 <- "=HYPERLINK(\"http://gnomad.broadinstitute.org/variant/"
    sGNOMAD2 <- "\",\"GNOMAD_link\")"
    variants <- variants %>%
        mutate(GNOMAD_POS = paste(Chrom, Pos, Ref, Alt, sep='-')) %>% 
        mutate(GNOMAD_Link = paste(sGNOMAD1, GNOMAD_POS, sGNOMAD2, sep = ''))

    # Columns 4,5: Ref,Alt

    # Column6 - Gene
    variants$Gene[variants$Gene == ""] <- NA
    
    # Column 6 - Zygosity, column 8 - Burden
    for(sample in samples){
        #DEBUG: 
        #gene = IL20RA
        sample <- samples[1]
        
        zygosity_column_name <- paste0("Zygosity.", sample)
        genotype_column_name <- paste0("gts.", sample)
        quo_test <- enquo("Ref")
        variants$Zygosity.WES_NA12878_202214 <- NULL
        variants <- variants %>% rename(!!genotype_column_name := test)
        
        variants <- variants %>% mutate(!!zygosity_column_name := !!genotype_column_name)
        (!!zygosity_column_name := str_length(!!quo_test))
        variants$test_zygosity <- NULL
                                        
        burden_column_name <- paste0("Burden.", sample)
        # calculating Burden using gene rather then Ensembl_gene_id - request from Matt
        t <- subset(variants, 
                    get(zygosity_column_name) == 'Hom' | get(zygosity_column_name) == 'Het',
                    select = c("Gene", zygosity_column_name))
        # count from plyr
        df_burden <- plyr::count(t, "Gene")    
        colnames(df_burden)[2] <- burden_column_name
        variants <- left_join(variants, df_burden, by = c("Gene"))
        variants[,burden_column_name][is.na(variants[,burden_column_name])] <- 0
        variants[,burden_column_name][is.na(variants$Gene)] <-0
    }
    
    # Column9 = gts
    
    # Column10 = Variation
    
    # Column11 =  Info
    variants <- add_column(variants, Info = rep("Info", nrow(variants)))
    
    for (i in 1:nrow(variants)){
        #debug: 
        i=1  
        v_id <- pull(variants, Variant_id)[i]
        gene <- variants[i,"Gene"]
        #for WES reports we need only coding impacts in the info field, for WGS we need all
        gene_impacts <- impacts %>%
            filter(variant_id == v_id) %>% 
            select(exon, hgvsc, hgvsp)
        if (coding) gene_impacts <- gene_impacts %>% filter(is_coding == 1)
        
        gene_impacts$gene <- rep(gene, nrow(gene_impacts))
        gene_impacts$exon[gene_impacts$exon==''] <- 'NA'
        gene_impacts <- gene_impacts %>% select(gene, exon, hgvsc, hgvsp)
        
        if (nrow(gene_impacts) > 0){
            v_impacts <- paste0(gene_impacts$gene, ":exon", gene_impacts$exon,
                                ":", gene_impacts$hgvsc, ":", gene_impacts$hgvsp)
            s_impacts <- paste(v_impacts, collapse = ",")
        }else s_impacts <- NA
      
        variants[i,"Info"] <- s_impacts
    }
    
    # Column12 - Refseq_change
    variants <- add_placeholder(variants, "Refseq_change", "NA")
    
    # Columns 13,14 - Depth, Quality

    # Column 15 - Alt_depth - from v.gt_alt_depths
    # when multiple callers used, AD is not set and fixed in merge_reports function
    for(sample in samples){
        new_name <- paste0("Alt_depths.", sample)
        setnames(variants, paste0("gt_alt_depths.", sample), new_name)
    }

    # Column 16 - Trio_coverage - fixed in merge_reports function
    variants <- add_placeholder(variants, "Trio_coverage", "")
    n_sample <- 1
    prefix <- ""
    
    #order gts column in the same way as in samples
    variants$gts <- ""
    for(sample in samples){
        column <- paste0("gt_depths.", sample)
        
        if (n_sample>1) prefix <- "/"
        
        variants$Trio_coverage <- with(variants, paste0(Trio_coverage, prefix, get(column)))
      
        column <- paste0("gts.", sample)
        
        if (n_sample>1) prefix <- ","
        
        variants$gts <- with(variants,paste0(gts, prefix, get(column)))
      
        n_sample <- n_sample+1
    }
    
    # Column17 = Ensembl_gene_id

    # Column18 = Gene_description
    gene_descriptions <- read_csv(paste0(default_tables_path, "/ensembl_genes_w_description.txt"))
    variants <- left_join(variants, gene_descriptions, by = c("Ensembl_gene_id" = "ensembl_gene_id"))
    
    # Column19 - Omim_gene_description
    omim_file_name <- paste0(default_tables_path,"/omim.txt")
    
    if (file.exists(omim_file_name)){
	    omim <- read.delim2(omim_file_name, stringsAsFactors=F)
	    variants <- merge(variants, omim, all.x=T)
	    variants$Omim_gene_description[is.na(variants$Omim_gene_description)] <- 0
    }
        
    # Column20 - Omim_inheritance 
    omim_inheritance_file_name <- paste0(default_tables_path,"/omim.inheritance.csv")        
    
    if (file.exists(omim_inheritance_file_name)){
	    omim_inheritance <- read.csv(omim_inheritance_file_name, stringsAsFactors = F)
	    variants <- merge(variants, omim_inheritance, all.x = T)
    }

    # Column 21 = Orphanet
    # previous name - orphanet.deduplicated.txt
    orphanet_file_name <- paste0(default_tables_path,"/orphanet.txt")
       
    if (file.exists(orphanet_file_name)){
	    orphanet <- read.delim(orphanet_file_name, stringsAsFactors = F)  
	    variants <- merge(variants, orphanet, all.x = T)
    
	    variants$Orphanet[is.na(variants$Orphanet)] = 0
    }
    
    # Column 22 - Clinvar
    
    # Column 23 - Ensembl_transcript_id
    
    # Column 24 - AA_position
    # changing separator from / to _ because otherwise excel converts it into date
    variants[,"AA_position"] <- with(variants, gsub("/", "_", AA_position), fixed = T)
    
    # Column 25 - Exon
    variants[,"Exon"] <- with(variants,gsub("/", "_", Exon), fixed = T)
    
    # Column 26 - Protein_domains
    
    # Column 27, 28 = Frequency_in_C4R, Seen_in_C4R_samples
    variants <- add_placeholder(variants, "Frequency_in_C4R", "Frequency_in_C4R")
    variants <- add_placeholder(variants, "Seen_in_C4R_samples", "Seen_in_C4R_samples")
    
    # Columns 29,30,31,32: HGMD
    for(hgmd_field in c("HGMD_id", "HGMD_gene", "HGMD_tag", "HGMD_ref")){
        variants <- add_placeholder(variants, hgmd_field, "NA")
    }

    # Column 33 - rsIds

    # population frequencies
    # Column34 = Gnomad_af
    # Column35 = Gnomad_af_popmax
    
    # Gnomad gene constraint scores
    # Column36 = Gnomad_oe_lof_score
    # Column37 = Gnomad_oe_mis_score
    gnomad_scores_file <- paste0(default_tables_path, "/gnomad_scores.csv")
    gnomad_scores <- read.csv(gnomad_scores_file, stringsAsFactors = F)
    variants <- merge(variants, gnomad_scores, all.x = T, all.y = F)

    # Column38 = Gnomad_ac
    # Column39 = Gnomad_hom
    for (field in c("Gnomad_ac","Gnomad_hom")){
        variants[,field] <- with(variants,gsub("-1", "0", get(field), fixed = T))
        variants[,field] <- with(variants,gsub("None", "0", get(field), fixed = T))
    }
    
    # Column41 - Conserved_in_20_mammals
    
    # pathogenicity scores
    # Column42 = sift
    # Column43 = polyphen
    # Column44 = cadd
    # Column45 = vest3
    for (i in 1:nrow(variants)){
        v_vest <- strsplit(variants[i,"Vest3_score"], ",", fixed = T)[[1]]
        variants[i, "Vest3_score"] <- max(v_vest)
    }
    
    # Column45 = revel
    
    # Column46 = Gerp
    
    # Column47 = Imprinting_status
    # Column48 = Imprinting_expressed_allele
    imprinting_file_name <- paste0(default_tables_path, "/imprinting.txt")
    imprinting <- read.delim(imprinting_file_name, stringsAsFactors = F)
    variants <- merge(variants, imprinting, all.x = T)
    
    # Column49 - pseudoautosomal
    pseudoautosomal_file_name <- paste0(default_tables_path, "/pseudoautosomal.txt")
    pseudoautosomal <- read.delim(pseudoautosomal_file_name, stringsAsFactors = F)
    variants <- merge(variants, pseudoautosomal, all.x = T)
    
    # Column50 - splicing
    variants <- add_placeholder(variants, "Splicing", "NA")
    if ("spliceregion" %in% colnames(impacts))
    {
        for (i in 1:nrow(variants)){
	          v_id <- variants[i,"Variant_id"]
	          splicing_impacts <- subset(impacts, variant_id == v_id,
	                                    select = c("maxentscan_diff","spliceregion"))
	          splicing_impacts <- subset(splicing_impacts, !is.na(maxentscan_diff))
	          splicing_impacts <- unique(splicing_impacts[order(splicing_impacts$maxentscan_diff),])
	          # capture the absolute difference - very weak site, or very strong site
	          # negative - strong alt, + weak alt.
	
	          s_splicing_field <- 0
	          
	          if (nrow(splicing_impacts) > 0){
	              strongest_alt_site <- head(splicing_impacts, n=1)
	              s_splicing_field <- strongest_alt_site$maxentscan_diff
	          }
	          
	          if (nrow(splicing_impacts) > 1){
	              weakest_alt_site <- tail(splicing_impacts, n=1)
	              s_splicing_field <- paste0(s_splicing_field, ";", weakest_alt_site$maxentscan_diff)
	          }
	          
	          variants[i,"Splicing"] <- s_splicing_field
        }
    }else print("VEP MaxEntScan annotation is missing")
    
    # Column 51: number of callers
    variants <- add_placeholder(variants, "Number_of_callers", "Number_of_callers")
    
    # Column 52: Old multiallelic
    variants$Old_multiallelic[variants$Old_multiallelic == "None"] <- "NA"
        
    # replace -1 with 0
    for (field in c("Trio_coverage", "Gnomad_af", "Gnomad_af_popmax")){
        variants[,field] <- with(variants, gsub("-1", "0", get(field), fixed = T))
        variants[,field] <- with(variants, gsub("None", "0", get(field), fixed = T))
    }

    for (field in c(paste0("Alt_depths.",samples))){
        variants[,field] <- with(variants, gsub("-1", NA, get(field), fixed = T))  
    }

    print(sort(colnames(variants)))
    select_and_write2(variants, samples, paste0(family, ".create_report"))
}

# writes in CSV format
select_and_write2 <- function(variants, samples, prefix)
{
    variants <- variants[c(c("Position", "UCSC_Link", "GNOMAD_Link", "Ref", "Alt"),
                          paste0("Zygosity.", samples),
                          c("Gene"),
                          paste0("Burden.", samples),
                          c("gts", "Variation", "Info", "Refseq_change", "Depth", "Quality"),
                          paste0("Alt_depths.", samples),
                          c("Trio_coverage", "Ensembl_gene_id", "Gene_description", "Omim_gene_description", "Omim_inheritance",
                            "Orphanet", "Clinvar",
                            "Frequency_in_C4R", "Seen_in_C4R_samples", "HGMD_id", "HGMD_gene", "HGMD_tag", "HGMD_ref",
                            "Gnomad_af_popmax", "Gnomad_af", "Gnomad_ac", "Gnomad_hom",
                            "Ensembl_transcript_id", "AA_position", "Exon", "Protein_domains", "rsIDs",
                            "Gnomad_oe_lof_score", "Gnomad_oe_mis_score", "Exac_pli_score", "Exac_prec_score", "Exac_pnull_score",
                            "Conserved_in_20_mammals", "Sift_score", "Polyphen_score", "Cadd_score", "Vest3_score", "Revel_score", "Gerp_score",
                            "Imprinting_status", "Imprinting_expressed_allele", "Pseudoautosomal", "Splicing",
                            "Number_of_callers", "Old_multiallelic"))]
  
    variants <- variants[order(variants$Position),]
    
    write.csv(variants, paste0(prefix,".csv"), row.names = F)
}

fix_column_name <- function(column_name){
    if(grepl("^[0-9]", column_name)){
        column_name <- paste0("X", column_name)
    }
    return(column_name)
}

# merges ensembl, gatk-haplotype reports
merge_reports <- function(family, samples){
    ensemble_file <- paste0(family, ".create_report.csv")
    ensemble <- read.csv(ensemble_file, stringsAsFactors = F)
    ensemble$superindex <- with(ensemble, paste(Position, Ref, Alt, sep = '-'))
    
    for (i in 1:nrow(ensemble)){
        v_impacts <- strsplit(ensemble[i,"Info"], "," , fixed = T)[[1]]
	    for (impact in v_impacts){
            if (grepl(":NM_", impact, fixed = T)){
                v_subimpacts <- strsplit(impact, ":", fixed=T)[[1]]
                ensemble[i,"Refseq_change"] <- paste0(v_subimpacts[3], ":", v_subimpacts[4], ":", v_subimpacts[6])
                break
            }
        }
    }
    
    ensemble_table_file <- paste0(family, ".table")
    if (file.exists(ensemble_table_file)){
        ensemble_table <- read.delim(ensemble_table_file, stringsAsFactors = F)
        ensemble_table$superindex <- with(ensemble_table, paste(paste0(CHROM,":",POS), REF, ALT, sep = '-'))
        ensemble_table[c("CHROM", "POS", "REF", "ALT")] <- NULL
        for (i in 1:nrow(ensemble_table)){
            if(!is.na(ensemble_table[i, "CALLERS"])){
        	      v_callers <- strsplit(ensemble_table[i, "CALLERS"],",")[[1]]
        	      ensemble_table[i, "Number_of_callers"] <- length(v_callers)
            }else ensemble_table[i,"Number_of_callers"] <- NA
        }
        ensemble_table["CALLERS"] <- NULL
        ensemble$Number_of_callers <- NULL
        #two variant callers called one genotype, two another - two genotypes, creates two records at the same site
        ensemble <- merge(ensemble, ensemble_table, by.x = "superindex", 
                          by.y = "superindex",all.x = T, all.y = F)
    }

    gatk_file <- paste0(family,"-gatk-haplotype-annotated-decomposed.table")
    if (file.exists(gatk_file)){
        gatk <- read.delim(gatk_file, stringsAsFactors = F)
        gatk$superindex <- with(gatk, paste(paste0(CHROM, ":", POS), REF, ALT, sep = '-'))
        gatk[c("CHROM","POS","REF","ALT")] <- NULL
    
        ensemble <- merge(ensemble, gatk, by.x = "superindex", by.y = "superindex", all.x = T, all.y = F)
    
        ensemble$Depth <- ensemble$DP
        n_sample <- 1
        prefix <- ""
        ensemble$Trio_coverage <- ""
    
        for(sample in samples){
            #R fixes numerical column names with X?
            #what if sample is not numerical
            column <- fix_column_name(sample)  
            column <- paste0(column,".DP")
        
	        #prefix changed to _ from / because otherwise excel converts the field into date
            if (n_sample > 1) prefix <- "_"
        
            ensemble$Trio_coverage <- with(ensemble, paste0(Trio_coverage, prefix, get(column)))
      
            column <- paste0("Alt_depths.", sample)
            column_gatk <- fix_column_name(sample)
            column_gatk <- paste0(column_gatk, ".AD")
        
            ensemble[,column] <- ensemble[,column_gatk]
      
            n_sample <- n_sample + 1
        }
    
        for (i in 1:nrow(ensemble)){
            for (sample in samples){
                field <- paste0("Alt_depths.", sample)
                #when combining reports from vcfs called elsewere there may be no AD field, just -1
                if (grepl(",", ensemble[i,field])){
            	    ensemble[i, field] <- strsplit(ensemble[i,field], ",", fixed = T)[[1]][2]
            	}
            }
        }
    
        for (sample in samples){
            ensemble[c("DP", paste0(fix_column_name(sample), ".DP"), 
                       paste0(fix_column_name(sample),".AD"))] <- NULL
        }
    }

    freebayes_file <- paste0(family,"-freebayes-annotated-decomposed.table")
    if(file.exists(freebayes_file)){
        freebayes <- read.delim(freebayes_file, stringsAsFactors = F)
        freebayes$superindex <- with(freebayes, paste(paste0(CHROM,":",POS), REF, ALT, sep = '-'))
        freebayes[c("CHROM","POS","REF","ALT")] <- NULL

        ensemble <- merge(ensemble, freebayes, by.x = "superindex", 
                          by.y = "superindex", all.x = T, all.y = F)
        for (i in 1:nrow(ensemble)){
            #if(grepl("NA",ensemble[i,"Trio_coverage"]))
            #wrong: a variant may be called by gatk with 10/10/NA,
            #and freebayes will destroy coverage info
            if (str_count(ensemble[i,"Trio_coverage"], "NA") == length(samples)){
                ensemble[i, "Depth"] <- ensemble[i, "DP"]
                for (sample in samples){
                    field_depth <- paste0("Alt_depths.", sample)
                    field_bayes <- paste0(fix_column_name(sample), ".AO")
                    #field_bayes = paste0(sample,".AO")
                    ensemble[i, field_depth] <- ensemble[i, field_bayes]
                }
                n_sample <- 1
                prefix <- ""
                ensemble[i, "Trio_coverage"] <- ""
            
                for(sample in samples){
                    column <- paste0(fix_column_name(sample),".DP")
                    if (n_sample > 1) prefix <- "_"
                        ensemble[i, "Trio_coverage"] <- paste(ensemble[i,"Trio_coverage"],
                                                              ensemble[i,column], sep = prefix)
              
                    n_sample <- n_sample+1
                }
          }
      }
      for (sample in samples){
          ensemble[c("DP", paste0(fix_column_name(sample),".DP"), 
                     paste0(fix_column_name(sample),".AO"))] <- NULL
      }
    }

    platypus_file <- paste0(family, "-platypus-annotated-decomposed.table")
    if(file.exists(platypus_file)){
        platypus <- read.delim(platypus_file, stringsAsFactors = F)
        platypus$superindex <- with(platypus, paste(paste0(CHROM,":",POS), REF, ALT, sep = '-'))
        platypus[c("CHROM", "POS", "REF", "ALT")] <- NULL
        ensemble <- merge(ensemble, platypus, by.x = "superindex", by.y = "superindex", 
                          all.x = T, all.y = F)
    
        for (i in 1:nrow(ensemble)){
          #if(grepl("NA",ensemble[i,"Trio_coverage"])) - wrong, may be 10/10/NA in gatk
          #if (ensemble[i,"Trio_coverage"]=="NA/NA/NA")
            if (str_count(ensemble[i,"Trio_coverage"],"NA") == length(samples)){
                ensemble[i,"Depth"] <- ensemble[i,"TC"]
                for (sample in samples){
                    field_depth <- paste0("Alt_depths.", sample)
                    field_bayes <- paste0(fix_column_name(sample), ".NV")
          
                    #sometimes freebayes has 10,10,10 for decomposed alleles
                    if (grepl(",", ensemble[i,field_bayes])){
                        ensemble[i,field_depth] <- strsplit(ensemble[i,field_bayes], ",", fixed = T)[[1]][1]
                    }
                }
                n_sample <- 1
                prefix <- ""
                ensemble[i, "Trio_coverage"] <- ""
        
                for(sample in samples){
                    column <- paste0(fix_column_name(sample), ".NR")
                    if (n_sample > 1) prefix <- "_"
                    #sometimes freebayes has 10,10,10 for decomposed alleles
                    if (grepl(",",ensemble[i,column])){
                        cov_value <- strsplit(ensemble[i,column], ",", fixed = T)[[1]][1]
                    }else cov_value <- ensemble[i,column]
                    
                    ensemble[i, "Trio_coverage"] <- paste(ensemble[i, "Trio_coverage"], cov_value, sep = prefix)
                    n_sample <- n_sample + 1
                }
          }
      }
    
      for (sample in samples){
          ensemble[c("TC", paste0(fix_column_name(sample), ".NV"), paste0(fix_column_name(sample),".NR"))] <-  NULL
      }
    }
    
    #don't use samtools file by default!
    samtools_file <- paste0(family,"-samtools-annotated-decomposed.table")
    if(file.exists(samtools_file)){
        samtools <- read.delim(samtools_file, stringsAsFactors = F)
        samtools$superindex <- with(samtools, paste(paste0(CHROM, ":", POS), REF, ALT, sep = '-'))
        samtools[c("CHROM", "POS", "REF", "ALT")] = NULL
        ensemble <- merge(ensemble, samtools, by.x = "superindex", 
                         by.y="superindex", all.x = T, all.y = F)
      
        for (i in 1:nrow(ensemble)){
            ensemble[i, "Depth"] = ensemble[i,"DP"]
            for (sample in samples){
                field_depth <- paste0("Alt_depths.", sample)
                field_samtools <- paste0(fix_column_name(sample), ".DP")
                ensemble[i, field_depth] <- ensemble[i, field_samtools]
            }
            ensemble[i, "Trio_coverage"] <- ""
        }
        for (sample in samples){
          ensemble[c("DP", paste0(fix_column_name(sample),".DP"))] <- NULL
          #samtools does not discriminate between insufficient coverage (cannot call) and no_call =reference
          field <- paste0("Zygosity.", sample)
          ensemble[,field] <- with(ensemble, gsub("Insufficient_coverage", 
                                   "-", get(field), fixed=T))
        }
    }
    
    ensemble[,"Trio_coverage"] <- with(ensemble,gsub("NA", "0", get("Trio_coverage"), fixed = T))  
   
    for (i in 1:nrow(ensemble)){
        if (is.na(ensemble[i, "Depth"])){
            l <- strsplit(ensemble[i, "Trio_coverage"],"_")[[1]]
            ensemble[i, "Depth"] <- sum(as.integer(l))
        }
        for (sample in samples){
            field_depth <- paste0("Alt_depths.", sample)
            if (is.na(ensemble[i, field_depth]))
                ensemble[i,field_depth] <- 0
        }
    }
    
    select_and_write2(ensemble, samples, paste0(family, ".merge_reports"))
}

annotate_w_care4rare <- function(family, samples){
    variants <- read.csv(paste0(family, ".merge_reports.csv"), stringsAsFactors = F)
  
    variants$superindex <- with(variants, paste(Position, Ref, Alt, sep='-'))
    
    if(exists("seen_in_c4r_counts")){
        variants <- merge(variants, seen_in_c4r_counts, by.x = "superindex", 
                          by.y = "Position.Ref.Alt", all.x = T)
        variants$Frequency_in_C4R <- variants$Frequency
        variants$Frequency <- NULL
    }
    
    variants$Frequency_in_C4R[is.na(variants$Frequency_in_C4R)] <- 0
    
    if(exists("seen_in_c4r_samples")){
        variants <- merge(variants,seen_in_c4r_samples,by.x = "superindex", 
                          by.y = "Position.Ref.Alt", all.x = T)
        variants$Seen_in_C4R_samples <- variants$Samples
    }
    
    variants$Seen_in_C4R_samples[is.na(variants$Seen_in_C4R_samples)] <- 0        
    
    if (exists("hgmd")){
        variants$HGMD_gene <- NULL
        variants$HGMD_id <- NULL
        variants$HGMD_ref <- NULL
        variants$HGMD_tag <- NULL
        variants <- merge(variants, hgmd, by.x = "superindex", 
                          by.y = "superindex", all.x = T, all.y = F)
        variants$HGMD_gene <- NULL
        
        hgmd.genes <- as.data.frame(unique(sort(hgmd$HGMD_gene)))
        hgmd.genes <- cbind(hgmd.genes, hgmd.genes)
        colnames(hgmd.genes) <- c("index", "HGMD_gene")
        variants <- merge(variants, hgmd.genes, by.x = "Gene", by.y = "index",
                          all.x = T, all.y = F)
    }
    
    select_and_write2(variants, samples, paste0(family, ".wes.", Sys.Date()))
}

load_tables <- function(debug = F){
    print(paste0("Debug:", debug))
    if (debug){
        seen_in_c4r_counts.txt <- "seen_in_c4r_counts.txt"
        seen_in_c4r_samples.txt <- "seen_in_c4r_samples.txt"
        hgmd.csv <- "hgmd.csv"
    }else{
        seen_in_c4r_counts.txt <- paste0(c4r_database_path, "/seen_in_c4r_counts.txt")    
        seen_in_c4r_samples.txt <- paste0(c4r_database_path, "/seen_in_c4r_samples.txt")
        hgmd.csv <- paste0(c4r_database_path, "/hgmd.csv")
    }
    
    if (file.exists(seen_in_c4r_counts.txt)){
        seen_in_c4r_counts <<- read_delim(seen_in_c4r_counts.txt, delim = "\t")
    }else print("No C4R counts found")
    
    if (file.exists(seen_in_c4r_samples.txt)){
        seen_in_c4r_samples <<- read_delim(seen_in_c4r_samples.txt, delim = "\t")
    }else print("No C4R samples found")
    
    if (file.exists(hgmd.csv)){
        hgmd <<- read_csv(hgmd.csv, 
                         col_names = c("chrom", "pos", "HGMD_id", "ref", "alt", "HGMD_gene",
                                       "HGMD_tag", "author", "allname","vol","page", "year", "pmid")) %>% 
                mutate(superindex = paste0(chrom, ':', pos, '-', ref, '-', alt),
                       HGMD_ref = paste(author, allname, vol, page, year, "PMID:", pmid, sep = ' ')) %>% 
                select(c(superindex, HGMD_id, HGMD_gene, HGMD_tag, HGMD_ref))
    }else print("No HGMD database")
}

# creates clinical report - more conservative filtering and less columns
clinical_report <- function(project, samples){
    report_file_name <- paste0(project,".wes.",Sys.Date(),".csv")
    
    full_report <- read_csv(report_file_name)
    #full_report <- mutate(full_report, max_alt = max(get(paste0("Alt_depths.", samples))))
    full_report$max_alt <- with(full_report, pmax(get(paste0("Alt_depths.", samples))))
    
    filtered_report <- subset(full_report, 
               Quality > 1000 & Gnomad_af_popmax < 0.005 & Frequency_in_C4R < 6 & max_alt >=20,
               select = c("Position", "GNOMAD_Link", "Ref", "Alt", "Gene", paste0("Zygosity.", samples), 
                          paste0("Burden.",samples),
                        "Variation", "Info", "Refseq_change", "Omim_gene_description", "Omim_inheritance",
                        "Orphanet", "Clinvar", "Frequency_in_C4R",
                        "Gnomad_af_popmax", "Gnomad_af", "Gnomad_ac", "Gnomad_hom",
                        "Sift_score", "Polyphen_score", "Cadd_score", "Vest3_score", "Revel_score",
                        "Imprinting_status", "Pseudoautosomal")
               )
    
    # recalculate burden using the filtered report
    for(sample in samples){
        zygosity_column_name <- paste0("Zygosity.", sample)
        burden_column_name <- paste0("Burden.", sample)
        t <- subset(filtered_report, 
                    get(zygosity_column_name) == 'Hom' | get(zygosity_column_name) == 'Het',
                    select = c("Gene", zygosity_column_name))
        # count is from plyr
        df_burden <- count(t, "Gene")    
        colnames(df_burden)[2] <- burden_column_name
        filtered_report[,burden_column_name] <- NULL
        filtered_report <- merge(filtered_report, df_burden, all.x = T)
        filtered_report[,burden_column_name][is.na(filtered_report[, burden_column_name])] <- 0
        filtered_report[,burden_column_name][is.na(filtered_report$Gene)] <-0
    }
    
    #order columns
    filtered_report <- filtered_report[c("Position", "GNOMAD_Link", "Ref", "Alt", "Gene", paste0("Zygosity.", samples), 
      paste0("Burden.", samples),
      "Variation", "Info", "Refseq_change", "Omim_gene_description", "Omim_inheritance",
      "Orphanet", "Clinvar", "Frequency_in_C4R",
      "Gnomad_af_popmax", "Gnomad_af", "Gnomad_ac", "Gnomad_hom",
      "Sift_score", "Polyphen_score", "Cadd_score", "Vest3_score", "Revel_score",
      "Imprinting_status", "Pseudoautosomal")]

    write.csv(filtered_report, paste0(project, ".wes.clinical.", Sys.Date(), ".csv"), row.names = F)
}

library(tidyverse)
library(data.table) #what for?
default_tables_path <- "~/cre/data"
c4r_database_path <- "/hpf/largeprojects/ccm_dccforge/dccforge/results/database"

args <- commandArgs(trailingOnly = T)
family <- args[1]

coding <- if(is.null(args[2])) T else F

debug <- F

setwd(family)

samples <- read_lines("samples.txt")

load_tables(debug)
create_report(family, samples)
merge_reports(family,samples)
annotate_w_care4rare(family,samples)
clinical_report(family,samples)

setwd("..")