# variant report generator - not supported, old version for gemini load, no Gnomad WGS frequencies
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

# output : family.ensemble.txt
create_report = function(family,samples)
{
    file=paste0(family,"-ensemble.db.txt")
  
    variants = get_variants_from_file(file)

    #Column1 - Position
    variants$Position=with(variants,paste(Chrom,Pos,sep=':'))
    
    #Column2 - UCSC link
    sUCSC1="=HYPERLINK(\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hgt.out3=10x&position="
    sUCSC2="\",\"UCSC_link\")"
    variants$UCSC_Link=with(variants,paste(sUCSC1,Position,sUCSC2,sep=''))

    # Column3 = GNOMAD_Link
    variants$GNOMAD_POS = with(variants,paste(Chrom,Pos,Ref,Alt,sep='-'))
    sGNOMAD1 = "=HYPERLINK(\"http://gnomad.broadinstitute.org/variant/"
    sGNOMAD2 = "\",\"GNOMAD_link\")"
    variants$GNOMAD_Link = with(variants,paste(sGNOMAD1,GNOMAD_POS,sGNOMAD2,sep=''))

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
    
    # Column 10 -  Info
    variants = add_placeholder(variants,"Info","Info")
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
      
        variants[i,"Info"] = s_impacts
    }
    
    # Column 11 - Protein_change_ensembl
    
    # Column 12 - Protein_change_refseq
    variants = add_placeholder(variants,"Protein_change_refseq","NA")
    
    # Columns 13,14 - Depth, Qual_depth

    # Column 15 - Alt_depth - from v.gt_alt_depth
    # when multiple callers used, AD is not set and fixed in merge_reports function
    for(sample in samples)
    {
        new_name = paste0("Alt_depths.",sample)
        setnames(variants, paste0("gt_alt_depths.",sample),new_name)
    }

    # Column 16 - Trio_coverage - fixed in merge_reports function
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
    
    # Column 17 - Ensembl_gene_id

    # Column 18 - Gene_description
    gene_descriptions = read.delim2(paste0(default_tables_path,"/ensembl_w_description.txt"), stringsAsFactors=F)
    variants = merge(variants,gene_descriptions,by.x = "Ensembl_gene_id",by.y = "ensembl_gene_id",all.x=T)
    
    # Column 19 - Omim_gene_description
    omim_file_name = paste0(default_tables_path,"/omim.txt")    

    if (file.exists(omim_file_name))
    {
	   omim = read.delim2(omim_file_name, stringsAsFactors=F)
	   variants = merge(variants,omim,all.x=T)
	   variants$Omim_gene_description[is.na(variants$Omim_gene_description)] = 0
    }
        
    # Column 20 - Omim_inheritance 
    #omim_inheritance_file_name = paste0(default_tables_path,"/omim_inheritance.txt")
          
    #if (file.exists(omim_inheritance_file_name))
    #{
	#   omim_inheritance = read.csv(omim_inheritance_file_name, sep="",stringsAsFactors = F)
	#   variants = merge(variants,omim_inheritance,all.x=T)
    #}
    
    # Column20 - Omim_inheritance 
    omim_inheritance_file_name = paste0(default_tables_path,"/omim.inheritance.csv")        
    
    if (file.exists(omim_inheritance_file_name))
    {
        omim_inheritance = read.csv(omim_inheritance_file_name,stringsAsFactors = F)
        variants = merge(variants,omim_inheritance,all.x=T)
    }

    # Column 21 - Orphanet
    # previous name - orphanet.deduplicated.txt
    orphanet_file_name = paste0(default_tables_path,"/orphanet.txt")       
    
    if (file.exists(orphanet_file_name))
    {
	   orphanet = read.delim(orphanet_file_name, stringsAsFactors=F)  
	   variants = merge(variants,orphanet,all.x=T)
    
	   variants$Orphanet[is.na(variants$Orphanet)] = 0
    }
    
    # Column 22 - Clinvar
    
    # Column 23 - Ensembl_transcript_id
    
    # Column 24 - AA_position
    # changing separator from / to _ because otherwise excel converts it into date
    variants[,"AA_position"] = with(variants,gsub("/","_",AA_position),fixed=T)
    
    # Column 25 - Exon
    variants[,"Exon"] = with(variants,gsub("/","_",Exon),fixed=T)
    
    # Column 26 - Pfam_domain
    
    # Column 27, 28 = Frequency_in_C4R, Seen_in_C4R_samples
    variants = add_placeholder(variants,"Frequency_in_C4R","Frequency_in_C4R")
    variants = add_placeholder(variants,"Seen_in_C4R_samples","Seen_in_C4R_samples")
    
    # Column 29 - rsIds

    # Columns 30-35 - population frequencies

    # Columns 36-37, Exac scores
    exac_scores_file = paste0(default_tables_path,"/exac_scores.txt")
    exac_scores = read.delim(exac_scores_file, stringsAsFactors=F)
    variants = merge(variants,exac_scores,all.x=T)

    # Column 39 - Gnomad_het
    
    # Column 40 - Exac het
    
    # Column 41 - Gnomad_hom_alt
    
    # Column 40 - Conserved in 29 mammals instead of phastcons
    #https://www.biostars.org/p/150152/

    # Column 41-42-43: sift,polyphen,cadd scores

    
    # Columns 44,45 - imprinting
    imprinting_file_name = paste0(default_tables_path,"/imprinting.txt")
    imprinting = read.delim(imprinting_file_name, stringsAsFactors=F)
    variants = merge(variants,imprinting,all.x=T)
    
    # Column 46 - pseudoautosomal
    pseudoautosomal_file_name = paste0(default_tables_path,"/pseudoautosomal.txt")
    pseudoautosomal = read.delim(pseudoautosomal_file_name, stringsAsFactors=F)
    variants = merge(variants,pseudoautosomal,all.x=T)
    
    # Column 47 - splicing
    variants = add_placeholder(variants,"Splicing","NA")
    
    #in older runs (before Nov2017) there are no splicing fields in the database
    # we report the max vep_maxentscan_diff (ref-alt) score for a gene over all isoforms
    # https://github.com/Ensembl/VEP_plugins/blob/master/MaxEntScan.pm
    # that means the lesser the score (more in the negative), the stronger is Alt splice signal
    # diff = REF  - ALT
    if ("vep_spliceregion" %in% colnames(impacts))
    {
        for (i in 1:nrow(variants))
        {
	          v_id = variants[i,"Variant_id"]
	          splicing_impacts = subset(impacts,variant_id==v_id,
	                             select=c("vep_maxentscan_diff","vep_spliceregion"))
	          splicing_impacts = subset(splicing_impacts, !is.na(vep_maxentscan_diff))
	          splicing_impacts = unique(splicing_impacts[order(splicing_impacts$vep_maxentscan_diff),])
	          # capture the absolute difference - very weak site, or very strong site
	          # negative - strong alt, + weak alt.
	
	          s_splicing_field=0
	          
	          if (nrow(splicing_impacts) > 0)
	          {
	              strongest_alt_site = head(splicing_impacts,n=1)
	            
	              s_splicing_field = strongest_alt_site$vep_maxentscan_diff
	          }
	          
	          if (nrow(splicing_impacts) > 1)
	          {
	              weakest_alt_site = tail(splicing_impacts,n=1)
	              s_splicing_field = paste0(s_splicing_field,";",weakest_alt_site$vep_maxentscan_diff)
	          }
	          
	          variants[i,"Splicing"] = s_splicing_field
        }
    }
    
    # Column 48: number of callers
    # Column 49: genotypes of individual callers
    variants = add_placeholder(variants,"Number_of_callers","Number_of_callers")
    #variants = add_placeholder(variants,"Genotypes_of_callers","Genotypes_of_callers")
    
    
    # replace -1 with 0
    for (field in c("EVS_maf_aa","EVS_maf_ea","EVS_maf_all","Maf_1000g","Gnomad_maf","Maf_all","Gnomad_het","Exac_het","Gnomad_hom_alt","Trio_coverage"))
    {
        variants[,field] = with(variants,gsub("-1","0",get(field),fixed=T))  
    }

    for (field in c(paste0("Alt_depths.",samples)))
    {
        variants[,field] = with(variants,gsub("-1",NA,get(field),fixed=T))  
    }

    print(sort(colnames(variants)))
    select_and_write2(variants,samples,paste0(family,".create_report"))
}

# column selection and order
select_and_write = function(variants,samples,prefix)
{
    variants = variants[c(c("Position","UCSC_Link","GNOMAD_Link","Ref","Alt"),
                        paste0("Zygosity.",samples),c("Gene"),
                        paste0("Burden.",samples),c("gts","Variation","Info","Protein_change_ensembl","Protein_change_refseq","Depth","Quality"),
                        paste0("Alt_depths.",samples),
                        c("Trio_coverage","Ensembl_gene_id","Gene_description","Omim_gene_description","Omim_inheritance",
                          "Orphanet", "Clinvar","Ensembl_transcript_id","AA_position","Exon","Pfam_domain",
                          "Frequency_in_C4R","Seen_in_C4R_samples","rsIDs","Maf_1000g","EVS_maf_aa","EVS_maf_ea","EVS_maf_all",
                          "Gnomad_maf","Maf_all", "Exac_pLi_score","Exac_missense_score","Gnomad_het","Exac_het","Gnomad_hom_alt",
                          "Conserved_in_29_mammals","Sift_score","Polyphen_score","Cadd_score",
                          "Imprinting_status","Imprinting_expressed_allele","Pseudoautosomal","Splicing",
                          "Number_of_callers"))]
  
    write.table(variants,paste0(prefix,".txt"),quote=F,sep = ";",row.names=F)  
}

# writes in CSV format
select_and_write2 = function(variants,samples,prefix)
{
    variants = variants[c(c("Position","UCSC_Link","GNOMAD_Link","Ref","Alt"),paste0("Zygosity.",samples),c("Gene"),
                        paste0("Burden.",samples),c("gts","Variation","Info","Protein_change_ensembl","Protein_change_refseq","Depth","Quality"),
                        paste0("Alt_depths.",samples),c("Trio_coverage","Ensembl_gene_id","Gene_description","Omim_gene_description","Omim_inheritance",
                                                        "Orphanet", "Clinvar","Ensembl_transcript_id","AA_position","Exon","Pfam_domain",
                                                        "Frequency_in_C4R","Seen_in_C4R_samples","rsIDs","Maf_1000g","EVS_maf_aa","EVS_maf_ea","EVS_maf_all",
                                                        "Gnomad_maf","Maf_all", "Exac_pLi_score","Exac_missense_score","Gnomad_het","Exac_het","Gnomad_hom_alt",
                                                        "Conserved_in_29_mammals","Sift_score","Polyphen_score","Cadd_score",
                                                        "Imprinting_status","Imprinting_expressed_allele","Pseudoautosomal","Splicing",
                                                        "Number_of_callers"))]
  
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
    ensemble_file = paste0(family,".create_report.csv")
    ensemble = read.csv(ensemble_file, stringsAsFactors=F)
    ensemble$superindex=with(ensemble,paste(Position,Ref,Alt,sep='-'))
    
    for (i in 1:nrow(ensemble))
    {
        v_impacts = strsplit(ensemble[i,"Info"],",",fixed=T)[[1]]
	      for (impact in v_impacts)
        {
            if (grepl(":NP_",impact,fixed = T))
            {
                v_subimpacts = strsplit(impact,":",fixed=T)[[1]]
                ensemble[i,"Protein_change_refseq"] = paste0(v_subimpacts[5],":",v_subimpacts[6])
                break
            }
        }
    }
    
    ensemble_table_file = paste0(family,".table")
    if (file.exists(ensemble_table_file))
    {
        ensemble_table = read.delim(ensemble_table_file,stringsAsFactors = F)
        ensemble_table$superindex=with(ensemble_table,paste(paste0("chr",CHROM,":",POS),REF,ALT,sep='-'))
        ensemble_table[c("CHROM","POS","REF","ALT")]=NULL
        for (i in 1:nrow(ensemble_table))
        {
            if(!is.na(ensemble_table[i,"CALLERS"]))
            {
        	      v_callers = strsplit(ensemble_table[i,"CALLERS"],",")[[1]]
        	      ensemble_table[i,"Number_of_callers"] = length(v_callers)
    	       }
    	      else
    	      {
    		        ensemble_table[i,"Number_of_callers"] = NA
    	      }
        }
        ensemble_table["CALLERS"]=NULL
        ensemble$Number_of_callers=NULL
        #two variant callers called one genotype, two another - two genotypes, creates two records at the same site
        ensemble = merge(ensemble,ensemble_table,by.x="superindex",by.y="superindex",all.x=T, all.y=F)
    }
    
    gatk_file = paste0(family,"-gatk-haplotype-annotated-decomposed.table")
    if (file.exists(gatk_file))
    {
        gatk = read.delim(gatk_file, stringsAsFactors=F)
        gatk$superindex=with(gatk,paste(paste0("chr",CHROM,":",POS),REF,ALT,sep='-'))
        gatk[c("CHROM","POS","REF","ALT")]=NULL
    
        ensemble = merge(ensemble,gatk,by.x = "superindex", by.y="superindex",all.x = T,all.y=F)
    
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
        
	    #prefix changed to _ from / because otherwise excel converts the field into date
            if (n_sample>1) prefix="_"
        
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
                #when combining reports from vcfs called elsewere there may be no AD field, just -1
                if (grepl(",",ensemble[i,field]))
                {
            	    ensemble[i,field]=strsplit(ensemble[i,field],",",fixed=T)[[1]][2]
            	}
            }
        }
    
        for (sample in samples)
        {
            ensemble[c("DP",paste0(fix_column_name(sample),".DP"),paste0(fix_column_name(sample),".AD"))]=NULL
        }
    }
    
    freebayes_file = paste0(family,"-freebayes-annotated-decomposed.table")
    if(file.exists(freebayes_file))
    {
        freebayes = read.delim(freebayes_file, stringsAsFactors=F)
        freebayes$superindex=with(freebayes,paste(paste0("chr",CHROM,":",POS),REF,ALT,sep='-'))
        freebayes[c("CHROM","POS","REF","ALT")]=NULL
        ensemble = merge(ensemble,freebayes,by.x = "superindex", by.y="superindex",all.x = T,all.y=F)
    
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
                    if (n_sample>1) prefix="_"
                        ensemble[i,"Trio_coverage"] = paste(ensemble[i,"Trio_coverage"],ensemble[i,column],sep = prefix)
              
                    n_sample = n_sample+1
                }
          }
      }
    
      for (sample in samples)
      {
          ensemble[c("DP",paste0(fix_column_name(sample),".DP"),paste0(fix_column_name(sample),".AO"))]=NULL
      }
    }
    
    platypus_file = paste0(family,"-platypus-annotated-decomposed.table")
    if(file.exists(platypus_file))
    {
        platypus = read.delim(platypus_file, stringsAsFactors=F)
        platypus$superindex=with(platypus,paste(paste0("chr",CHROM,":",POS),REF,ALT,sep='-'))
        platypus[c("CHROM","POS","REF","ALT")]=NULL
        ensemble = merge(ensemble,platypus,by.x = "superindex", by.y="superindex", all.x = T, all.y = F)
    
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
                    if (n_sample>1) prefix="_"
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
    }
    
    #don't use samtools file by default!
    samtools_file = paste0(family,"-samtools-annotated-decomposed.table")
    if(file.exists(samtools_file))
    {
        samtools = read.delim(samtools_file, stringsAsFactors=F)
        samtools$superindex=with(samtools,paste(paste0("chr",CHROM,":",POS),REF,ALT,sep='-'))
        samtools[c("CHROM","POS","REF","ALT")]=NULL
        ensemble = merge(ensemble,samtools,by.x = "superindex", by.y="superindex",all.x = T, all.y=F)
      
        for (i in 1:nrow(ensemble))
        {
            ensemble[i,"Depth"] = ensemble[i,"DP"]
            for (sample in samples)
            {
                field_depth = paste0("Alt_depths.",sample)
                field_samtools = paste0(fix_column_name(sample),".DP")
                ensemble[i,field_depth] = ensemble[i,field_samtools]
            }
            ensemble[i,"Trio_coverage"]=""
        }
        for (sample in samples)
        {
          ensemble[c("DP",paste0(fix_column_name(sample),".DP"))]=NULL
          #samtools does not discriminate between insufficient coverage (cannot call) and no_call =reference
          field=paste0("Zygosity.",sample)
          ensemble[,field] = with(ensemble,gsub("Insufficient_coverage","-",get(field),fixed=T)) 
        }
    }
      
    
    ensemble[,"Trio_coverage"] = with(ensemble,gsub("NA","0",get("Trio_coverage"),fixed=T))  
   
    for (i in 1:nrow(ensemble))
    {
        if (is.na(ensemble[i,"Depth"]))
        {
            l=strsplit(ensemble[i,"Trio_coverage"],"_")[[1]]
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
        variants = merge(variants,seen_in_c4r_counts,by.x = "superindex", by.y="Position.Ref.Alt",all.x = T)
        variants$Frequency_in_C4R = variants$Frequency
        variants$Frequency=NULL
    }
    
    variants$Frequency_in_C4R[is.na(variants$Frequency_in_C4R)] = 0        
    
    if(exists("seen_in_c4r_samples"))
    {
        variants = merge(variants,seen_in_c4r_samples,by.x = "superindex", by.y="Position.Ref.Alt",all.x = T)
        variants$Seen_in_C4R_samples=variants$Samples
    }
    
    variants$Seen_in_C4R_samples[is.na(variants$Seen_in_C4R_samples)] = 0        
    
    if (exists("hgmd"))
    {
        variants = merge(variants,hgmd,by.x="superindex",by.y="superindex",all.x=T,all.y=F)
    }
    
    select_and_write2(variants,samples,paste0(family,".wes.",Sys.Date()))
}

library(stringr)
library(data.table)
library(plyr)

default_tables_path="~/cre/data"
c4r_database_path = "/hpf/largeprojects/ccm_dccforge/dccforge/results/database"

#load c4r information
seen_in_c4r_counts.txt = paste0(c4r_database_path,"/seen_in_c4r_counts.txt.chr")
if (file.exists(seen_in_c4r_counts.txt))
{
    seen_in_c4r_counts = read.delim(seen_in_c4r_counts.txt, stringsAsFactors=F)
}
seen_in_c4r_samples.txt = paste0(c4r_database_path,"/seen_in_c4r_samples.txt.chr")
if (file.exists(seen_in_c4r_samples.txt))
{
    seen_in_c4r_samples = read.delim(seen_in_c4r_samples.txt, stringsAsFactors=F)
}


# R substitutes "-" with "." in sample names in columns so fix this in samples.txt
# sample names starting with letters should be prefixed by X in *.table
# for correct processing. most of them start with numbers, and R adds X automatically

args = commandArgs(trailingOnly = T)
family = args[1]
debug = FALSE

setwd(family)

samples = unlist(read.table("samples.txt", stringsAsFactors=F))
samples = gsub("-",".",samples)
    
create_report(family,samples)
merge_reports(family,samples)
annotate_w_care4rare(family,samples)

setwd("..")

