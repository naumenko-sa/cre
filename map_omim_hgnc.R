# 01. Map HGNC Gene Names (Approved/Alternative/Previous) to OMIM Phenotypic Information and Inheritance
# Author: Delvin So
# 11/14/2019 - now reports multiple inheritances rather than just the first one
# 07/24/2020 - fixed so the final output is cast to wide, with each gene corresponding to a single string of omim phenotype and inheritances
# Parsing failures are from the comments at the end of the mim files
# Usage: Rscript map_omim_hgnc.R hgnc_website.txt mimTitles.txt genemap2.txt mim2gene.txt OMIM

suppressPackageStartupMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE) # returns argumentn after --args
options(warn = -1) # suppressing separate warnings
print("Warnings have been turned off to suppress tidyr::separate() warnings!")


hgnc_file_path <- args[1]
mim_titles_file <- args[2]
gene_map_file <- args[3]
mim2gene_file <- args[4]
output_prefix <- args[5]

print(hgnc_file_path)
print(mim_titles_file)
print(gene_map_file)
print(mim2gene_file)
print(output_prefix)
#mim_titles_file <- "mimTitles.txt"
#mim2gene_file <- "mim2gene.txt"
#gene_map_file <- "genemap2.txt"
# morbid_map_file <- "morbidmap.txt"
#hgnc_file_path <- "hgnc_website.txt"
#output_prefix <- Sys.Date()


# ---- read in HGNC and clean -----

hgnc_file <- read.csv(hgnc_file_path, 
                      sep = "\t",
                      stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rename_all(. %>% 
               str_to_lower() %>% 
               str_replace_all( "\\.", "_")) %>%
  rename(ensembl_id = "ensembl_gene_id", ncbi_id = "ncbi_gene_id", synonyms = "alias_symbols") %>%
  mutate(ncbi_id = as.character(ncbi_id), ensembl_id = as.character(ensembl_id))


# ---- read in OMIM phenotypic information and clean ----

# mimTitles.txt  ----

# Asterisk (*)  Gene
# Plus (+)  Gene and phenotype, combined
# Number Sign (#)  Phenotype, molecular basis known
# Percent (%)  Phenotype or locus, molecular basis unknown
# NULL (<null>)  Other, mainly phenotypes with suspected mendelian basis
# Caret (^)  Entry has been removed from the database or moved to another entry


titles_legend <- tibble(prefix = c("Asterisk", "Plus", "Number Sign", "Percent", "NULL", "Caret"),
                        symbol = c("*", "+", "#", "%", "NULL", "^"),
                        desc = c("gene", "gene + phenotype", "phenotype", "phenotype | locus", "Other - phenotype w/ suspected mendelian basis",
                                 "removed"))

mim_titles <- read_delim(mim_titles_file, delim = "\t", skip = 2,
                         col_names = c("prefix", "mim_number", "preferred_title", "alternative_titles", "included_titles")) %>%
  filter(!str_detect(prefix, "#")) %>%
  left_join(titles_legend, by = c("prefix"))


# mim2gene.txt ----
# A tab-delimited file linking MIM numbers with NCBI Gene IDs, Ensembl Gene IDs, and HGNC Approved Gene Symbols.
mim_gene <- read_delim(mim2gene_file, delim = "\t", skip = 4,
                       col_name = c("mim_number", "mim_entry_type", "ncbi_id", "hgnc_id", "ensembl_id"))  %>%
  # .[-c(1:5), ] %>%
  # set_names(c("mim_number", "mim_entry_type", "ncbi_id", "hgnc_id", "ensembl_id")) %>%
  mutate_at(vars(c("ncbi_id", "hgnc_id", "ensembl_id")), ~ ifelse(. == "", NA, .)) %>%
  as_tibble()

# morbidmap.txt (currently not used) ----

# A tab-delimited file of OMIM's Synopsis of the Human Gene Map (same as genemap2.txt above) sorted alphabetically by disorder.
# morbid_map <- read_delim(morbid_map_file, delim = "\t", skip = 3) %>%
#   as_tibble() %>%
#   set_names(c("phenotype", "gene_symbols", "mim_number", "cyto_loc")) %>%
#   filter(!str_detect(phenotype, "^#"))  
# 
# genemap2.txt ----
# A tab-delimited file containing OMIM's Synopsis of the Human Gene Map including additional information such as genomic coordinates and inheritance.

mim_map <- read_delim(gene_map_file, delim = "\t", skip = 3,
                      col_names = c("chromosome", "genomic_position_start", "genomic_position_end", "cyto_loc", "computed_cyto_loc",
                                    "mim_number", "gene_symbols", "gene_name", "approved_symbol", "entrez_gene_id", "ensembl_id", "comments",
                                    "phenotypes", "mouse_gene_symbol_id")) %>%
  # removes  metadata at bottom of the file
  filter(!str_detect(chromosome, "^#"))  


#----  start cleaning here ----

# ---- combine mim tables so we have entry type, phenotype, mim number, etc  -----


# join mim tables together to obtain phenotypes, gene symbols, names
mim_info <- mim_map %>%
  # left join to mim_gene b/c it has hgnc_id
  left_join(mim_gene, by = c("mim_number", "ensembl_id")) %>%
  left_join(mim_titles) %>%
  filter(str_detect(desc, "gene")) %>% # only retain genes, ~ 16187 + 37 / 16993
  mutate_all(~ ifelse(is.na(.), "", .)) %>%
  select(mim_number, gene_symbols, approved_symbol, ensembl_id, phenotypes) %>%
  mutate(complete = ifelse(!is.na(approved_symbol) & !is.na(ensembl_id), 1, 0))

# mim_info

# total_omim_genes <- str_detect(mim_info$desc, "gene")
# 
# print(paste0("Total Genes, Genes + Phenotype: ", total_omim_genes)) # debug, 15733

# 'flatten' the hgnc approved symbols/synonyms/previous symbols so each gene symbol is mapped to its hgnc_id
  
  # separate hgnc synonyms and convert into a stacked format
  hgnc_synonyms <- hgnc_file %>%
    select(hgnc_id, approved_symbol, synonyms, previous_symbols, ensembl_id, ncbi_id, status ) %>%
    # separate synonyms into a row for each synonym
    separate(synonyms, sep = ", ", into = paste0("syns_", c(1:15))) %>%
    # print(n = 100) %>%
    gather(key = "syns_name", "syns", contains("syns_"))  %>%
    filter(syns != "") %>% select(-(syns_name)) %>%
    select(hgnc_id, synonyms = syns)
  
  # separate previous hgnc previous symbols into a stacked format
  hgnc_prev <- hgnc_file %>% 
    select(hgnc_id, approved_symbol, synonyms, previous_symbols, ensembl_id, ncbi_id, status ) %>%
    # separate previous symbols into a row for each symbol
    separate(previous_symbols, sep = ", ", into = paste0("prev_", c(1:15)), extra = "drop") %>%
    gather(key = "prev_name", "prev", contains("prev_"))   %>%
    filter(prev != "") %>% select(-(prev_name)) %>%
    select(hgnc_id, previous_name = prev)

  # subset the hgnc approved symbols, w/ hgnc_ids
  hgnc_approv <- hgnc_file %>%
    select(hgnc_id, approved_symbol)
  
  # putting it all together
  all_hgnc_syns <- bind_rows(hgnc_approv %>%
                               rename(gene_name = approved_symbol),
            hgnc_synonyms %>%
              rename(gene_name = synonyms),
            hgnc_prev %>% 
              rename(gene_name = previous_name)) 
  
  # remove the (1) NAs
  all_hgnc_syns <- all_hgnc_syns %>% filter(!is.na(gene_name))

  print(paste0("Total # of HGNC Gene IDs: ",  length(unique(all_hgnc_syns$hgnc_id)))) # debug, 46612

 # ---- join mim phenotypes, joined to hgnc on approved symbols -----
# end result is each mim number is mapped to an hgnc_id using 'approved_symbol' (99.3% match rate)

mim_join_hgnc <- mim_info %>%
  left_join(hgnc_approv, by = c("approved_symbol")) 
  
#nrow(mim_info) < nrow(mim_join_hgnc) b/c symbols map to multiple hgnc id
  # mim_join_hgnc %>% count(approved_symbol, sort = TRUE)
  
# sanity check
# mim_join_hgnc %>% count(!is.na(hgnc_id)) %>% mutate(perc = prop.table(n) * 100)
# 99.3 % 
# mim_join_hgnc %>% filter(is.na(hgnc_id))

# ----  what about the 0.7% of OMIM that did not map to an HGNC ID? ----

# mim_join_hgnc %>% filter(is.na(hgnc_id)) %>% write_csv("output/mim_genes_unmapped_with_info.csv")

# mim_info %>% 
#   filter(mim_number %in% c(mim_join_hgnc  %>% filter(is.na(hgnc_id)) %>% 
#   select(mim_number) %>% pull(mim_number))) 

# these do not have any gene symbols..

# flatten the gene symbols..
mim_unmapped_hgnc <- mim_join_hgnc %>%
  filter(is.na(hgnc_id)) %>%
  separate(gene_symbols, sep = ", ", into = paste0("gene_", c(1:10)), extra = "drop") %>%
  # naniar::miss_var_summary() # debug
  # print(n = 100) %>%
  gather(key = "gene_name", "gene_symbol", contains("gene_")) %>% 
  select(mim_number, gene_symbol) 

# mim_unmapped_hgnc

# try to map the now flattened omim gene symbols onto hgnc
mim_unmapped_hgnc2 <- mim_unmapped_hgnc  %>%
  left_join(all_hgnc_syns, by = c("gene_symbol" = "gene_name")) 

# how many ended up mapping?
# mim_unmapped_hgnc2 %>%
#   group_by(mim_number) %>%
#   arrange(hgnc_id) %>%
#   print
# ---- some mim_ids map to multiple HGNC ids!!!!!!!!!! ----
# eg. BCL5 is found in HGNC 1000 and 1001
# 
# all_hgnc_syns %>% filter(gene_name == "BCL5")
# mim_info %>% filter(mim_number == 151441)

# mim_unmapped_hgnc2 %>% distinct(mim_number, hgnc_id) %>% na.omit() %>% arrange(mim_number) %>% print(n = 100)

# take the first hgnc id for each mim number because multiple gene symbols exist for each mim and thus can map to the same hgnc id more than once
mim_mapped_unmapped_hgnc <- mim_unmapped_hgnc2 %>%
  group_by(mim_number) %>%
  arrange(hgnc_id) %>%
  filter(row_number() == 1) %>% ungroup()

mim_unmapped_hgnc2  %>% distinct(mim_number, hgnc_id) %>% na.omit

# #how many still remain unmapped?
# mim_mapped_unmapped_hgnc %>%
#   count(!is.na(hgnc_id)) %>%
#   mutate(prop = prop.table(n))
# #68/112 or 60%

# the gene symbol is intermediary, the end goal is to map the mim number to an hgnc id
# mim_unmapped_hgnc2 %>%
#   distinct(mim_number, hgnc_id) %>%
#   na.omit() 

# join the unmapped IDs 
mim_join_hgnc2 <- mim_join_hgnc %>%
  # filter(!mim_number %in% c(   mim_unmapped_hgnc2 %>% 
  #                                                         distinct(mim_number) %>% pull(mim_number))) %>%
  left_join(mim_unmapped_hgnc2 %>%
              distinct(mim_number, hgnc_id) %>% 
              na.omit() ,
            by = c("mim_number")) %>%
  mutate(hgnc_id = coalesce(hgnc_id.x, hgnc_id.y))


# mim_join_hgnc2  %>% print(n = 100)
# # sanity check, the total between x and y should be the total in hgnc_id
# mim_join_hgnc2 %>%
#   summarize_at(vars(hgnc_id.x, hgnc_id.y, hgnc_id), ~ sum(!is.na(.)))

# A tibble: 1 x 3
# hgnc_id.x hgnc_id.y hgnc_id
# <int>     <int>   <int>
#   1     15633        83   15716
# looks good

# remove the old hgnc id columns
mim_join_hgnc2 <- mim_join_hgnc2 %>% select(-c(hgnc_id.x, hgnc_id.y))

# remaining unmapped

print("Total and Percentage of OMIM Genes Mapped to HGNC")
# print("")
totals <- mim_join_hgnc2 %>% 
  # distinct() %>%
  mutate(mapped = ifelse(is.na(hgnc_id), "Unmapped", "Mapped")) %>% 
  count(mapped) %>% mutate(prop = round(prop.table(n) * 100, 2)) %>%
  set_names(c("# OMIM Mapped to HGNC", "Count"), "Percentage") %>%
  data.frame()
totals[2, ] <- c("Totals", sum(totals$Count), sum(totals$Percentage))
print(totals)

  # spread(`# OMIM Mapped to HGNC`, Count)  # debug, 15716/44
print("Total # of genes above may not match with total # of genes in OMIM because an individual gene symbol # may have mapped to multiple mim_numbers. eg. approved_symbol = XAGE1B")

# ---- save the unmapped gene symols ----

# unmapped_file_name <- paste0("output/", output_prefix, "_remaining_unmapped_omim_with_info_", Sys.Date(), ".tsv")
unmapped_file_name <- paste0(output_prefix, "_remaining_unmapped_omim_with_info_", Sys.Date(), ".tsv")

mim_info %>% 
  # remove mim numbers that were NA (ie. unmapped to HGNC )
  filter(mim_number %in% c(mim_mapped_unmapped_hgnc %>% 
                             filter(is.na(hgnc_id)) %>% pull(mim_number) )) %>% # should be 44 (out of 112)
  # print(n = 50) %>%
  write_tsv(unmapped_file_name)

print(paste0("List of OMIM genes unable to be mapped to HGNC saved to: ", unmapped_file_name))

# ----- now, join the mim ids,, which have been mapped to hgnc ids, with all possible gene names on hgnc_id  -----

mim_join_hgnc_res <- mim_join_hgnc2 %>%
  # select(-c(ensembl_id, complete)) %>%
  
  # convert NA hgnc id to "" so it doesn't join on the NAs in all_hgnc_syns
  mutate(hgnc_id = ifelse(is.na(hgnc_id), "", hgnc_id)) %>%
  # the join
  left_join(all_hgnc_syns, by = "hgnc_id")  %>%
  # add the unmapped gene names back
  # filter(!mim_number %in% c(   mim_unmapped_hgnc2 %>% 
  #                                distinct(mim_number) %>% pull(mim_number))) %>%
  # bind_rows(mim_unmapped_hgnc2) %>%
  # re-convert all ""'s to NAs
  mutate_all(~ ifelse(. == "", NA, .)) %>%
  arrange(phenotypes) %>%
  # count(approved_symbol == gene_name) %>% # sanity check, looks good
  select(-c(gene_symbols, approved_symbol))

# mim_join_hgnc_res
# ----- extract phenotype and inheritance from 'phenotypes' column ----
# new code accounting for all inheritance values
hgnc_join_omim_phenos <- mim_join_hgnc_res %>% 
  # filter(!is.na(phenotypes)) %>% 
  # select(hgnc_id, approved_name, mim_number, phenotypes, ensembl_id) %>%
  mutate(phenotypes = str_to_lower(phenotypes)) %>%
  # this grabs everything after any occurence of and including 'Autosomal|X-Linked|Y-Linked' up until a semi colon, comma or nothing
  mutate(mim_inheritance = str_extract_all(phenotypes, "(autosomal|x-linked|y-linked).*?(\\;|\\,|$)")) %>% 
  unnest(mim_inheritance) %>%
  group_by(mim_number) %>%
  # collapse 
  summarize(mim_inheritance = paste0(str_to_lower(mim_inheritance), collapse = ",")) %>%
  ungroup()  %>%
  distinct()


hgnc_join_omim_phenos <- mim_join_hgnc_res %>%
  left_join(hgnc_join_omim_phenos %>% select(mim_number, mim_inheritance), by = "mim_number") %>%
  rename(mim_inheritance_extract = "mim_inheritance") %>% 
  # replace the first occurence of each inheritance with abbreviation, the remaining ones can be removed so we retain only unique occurences
  mutate(mim_inheritance = str_replace(mim_inheritance_extract, "autosomal recessive", "AR"),
         mim_inheritance = str_replace(mim_inheritance, "autosomal dominant", "AD"),
         mim_inheritance = str_replace(mim_inheritance, "x(\\-|\\s)linked recessive", "XL-R"),
         mim_inheritance = str_replace(mim_inheritance, "x(\\-|\\s)linked dominant", "XL-D"),
         mim_inheritance = str_replace(mim_inheritance, "x(\\-|\\s)linked", "XL"),
         mim_inheritance = str_replace(mim_inheritance, "y(\\-|\\s)linked recessive", "YL-R"),
         mim_inheritance = str_replace(mim_inheritance, "y(\\-|\\s)linked dominant", "YL-D"),
         mim_inheritance = str_replace(mim_inheritance, "y(\\-|\\s)linked", "YL"),
         mim_inheritance = str_replace_all(mim_inheritance, "[0-9]*", ""),  # removes numbers and lowercase letters (chromosomal?)
         mim_inheritance = str_replace_all(mim_inheritance, "NA(,*)", ""),
         mim_inheritance = str_replace_all(mim_inheritance, "[a-z]*", "")) %>%
  # remove anything that was not abbreviated since they must be duplicates as well as symbols
  mutate(mim_inheritance = str_replace_all(mim_inheritance, "(autosomal recessive)|(autosomal dominant)|x(\\-|\\s)linked|y(\\-|\\s)linked|;|-|\\(|\\)|\\{|\\}|\\\\|\\/", ""),
                mim_inheritance = str_replace_all(mim_inheritance, ",+", ","), # sub all instances of multipe commas with 1 comma
                mim_inheritance = str_replace_all(mim_inheritance, "\\s+", ""), # sub all white space with 1 
                mim_inheritance = str_replace_all(mim_inheritance, "\\B,\\s*", "")) # remove trailing commas bounded


# ----- new additions from 07/2020 -----
# cast so each row is a gene with omim description by casting to wide and concatenating together
final <- hgnc_join_omim_phenos  %>%  
  #.[1:1000,] %>% 
  select(omim_inheritance = mim_inheritance, omim_phenotype = phenotypes, gene_name)  %>% 
  mutate(omim_phenotype = str_to_lower(omim_phenotype)) %>% 
  group_by(gene_name) %>% 
  #filter(gene_name == 'RP1') %>%
  mutate(id = row_number()) %>% 
  gather(key = 'variable', value = 'value', c(omim_inheritance, omim_phenotype)) %>% 
  mutate(var_id = paste0(variable, id)) %>%
  select(-variable, -id) %>% 
  spread(var_id, value) %>%
  tidyr::unite_("omim_phenotype", sep = ",", names(.)[str_detect(names(.), 'phenotype')]) %>%
  tidyr::unite_("omim_inheritance", sep = ",", names(.)[str_detect(names(.), 'inheritance')]) %>% 
  mutate_at(vars(omim_inheritance, omim_phenotype), ~  str_replace_all(., "NA", "") %>% 
              str_replace_all("autosomal|recessive|autosomal|dominant|x(\\-|\\s)linked|y(\\-|\\s)linked", "") %>% 
              str_replace_all("(,|\\;|\\?)+", ",") %>% 
              str_replace_all("\\(|\\)|\\{|\\}|\\\\|\\/", "") %>% 
              str_replace_all("[0-9]*", "") %>% 
              str_replace_all("\\s+", " ") %>%
              str_replace_all("\\B,\\s*", "") %>%
              str_replace_all("(,$)|(,\\s$)", "")) %>%
  mutate_at(vars(omim_inheritance, omim_phenotype), ~ ifelse(.x == "", NA, .x)) %>%
  arrange(omim_inheritance)
# tidyr::unite_("Omim", sep = ";", names(.)[-1]) %>% 
  # mutate(Omim = str_replace_all(Omim, "NA", ""),
  #        Omim = str_replace_all(Omim, "autosomal|recessive|autosomal|dominant|x(\\-|\\s)linked|y(\\-|\\s)linked", ""),
  #        Omim = str_replace_all(Omim, "(,|\\;|\\?)+", ","),
  #        Omim = str_replace_all(Omim, "\\(|\\)|\\{|\\}|\\\\|\\/", ""),
  #        Omim = ifelse(Omim == ",", NA, Omim),
  #        Omim = str_replace_all(Omim, "[0-9]*", ""),  # removes numbers and lowercase letters (chromosomal?))
  #        Omim = str_replace_all(Omim, "\\s+", " "), # sub all white space with 1 
  #        Omim = str_replace_all(Omim, "\\B,\\s*", ""))

stopifnot(length(unique(final$gene_name)) == nrow(final))
final
# ---- save to output file -----
hgnc_omim_output_file_name <- paste0(output_prefix, "_hgnc_join_omim_phenos_", Sys.Date(), ".tsv")

write_tsv(final, hgnc_omim_output_file_name)
print(paste0("Output saved to: ", hgnc_omim_output_file_name))


  


  
  
