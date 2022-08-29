# author: "Elsa Leitao"

library(readr)
library(dplyr)
library(readxl)
library(stringr)
library(tidyr)


### Genes with approved symbols from HGNC 
############################################################################################################################

# The genes on all chromosomes were downloaded from HGNC ("data/HGNC_symbols_allChr.txt"). 
# Search done in HGNC (https://www.genenames.org/download/custom/) 

## Curated by the HGNC
# HGNC ID
# Approved symbol
# Approved name
# Status
# Locus type
# Locus group
# Previous symbols
# Alias symbols
# Chromosome
# Ensembl ID

## Downloaded from external sources
# NCBI Gene ID (supplied by NCBI)
# Uniprot ID (supplied by Uniprot)
# Ensembl ID (supplied by Ensembl)
# Mouse genome database ID (supplied by MGI)
# Rat genome database ID (supplied by RGD)

## Select status
# Approved
# Entry and symbol withdrawn


HGNC_symbols_allChr <- read_delim("data/HGNC_symbols_allChr.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  rowwise() %>%
  mutate(chr = str_split(Chromosome, "p|q|[:blank:]") %>% unlist() %>% .[1]) %>%
  mutate(chr = factor(chr, levels = c(as.character(1:22), "X", "Y")))

HGNC_symbols_allChr_approved <- HGNC_symbols_allChr %>%
  filter(all(Status == "Approved", `Locus group` == "protein-coding gene", !grepl("not on reference assembly", Chromosome)))


### OMIM annotations   
############################################################################################################################

# genemap2.txt file from OMIM (https://www.omim.org/downloads) 

names_col = read_tsv("data/OMIM_genemap2.txt",
                     col_names = FALSE,
                     skip = 3,
                     n_max = 1,
                     col_types = cols())

omim_genes_hg38 = read_tsv("data/OMIM_genemap2.txt",
                           col_names = FALSE,
                           comment = "#",
                           col_types = cols()) 

names(omim_genes_hg38) = names_col

omim_genes_hg38_annot = omim_genes_hg38 %>%
  rename(Chromosome = "# Chromosome") %>%
  rowwise() %>%
  mutate(Chromosome = str_replace(Chromosome, "chr", "")) %>%
  mutate(Chromosome = factor(Chromosome, levels = c(as.character(1:22), "X", "Y"))) %>%
  mutate(morbid = case_when(!is.na(`Approved Symbol`) & !is.na(Phenotypes) ~ T,
                            T ~ F))


## OMIM morbid genes (subset of genes that have "Approved Symbol" & "Phenotype")

omim_morbid_all = omim_genes_hg38_annot %>% 
  filter(morbid == T)


## Separation of multiple phenotypes into different rows
## and split into 3 columns the phenotype name, MIM number and mode of inheritance.  

omim_morbid_split = omim_morbid_all %>%
  separate_rows(Phenotypes, sep = "; ") %>%
  rowwise() %>%
  mutate(phen = str_split(Phenotypes, ", [:digit:]{6}") %>% unlist %>% .[1],
         phen_MIM_number = str_extract(Phenotypes, "[:digit:]{6}"),
         inheritance = str_split(Phenotypes, "[(][1-4][)], |[(][1-4][)]") %>% unlist %>% last()) %>%
  mutate(inheritance = case_when(inheritance == "" ~ NA_character_,
                                 TRUE ~ inheritance))


## Download Clinical Synopsis data for OMIM phenotypic entries

# OMIM Advanced Search for entries with (restricting to "MIM Number Prefix: # phenotype description, molecular basis known"):  
# - Clinical Synopsis  
# - Clinical Synopsis "neurologic" field (+cs_neurologic_exists:true)  
# 
# OMIM Advanced Search for entries with terms in the Clinical Synopsis "neurologic" field related to:  
# - intelectual disability  
# - impaired motor development  
# - impaired language development  
# - seizures  
# - spasticity  
# - ataxia  
# 
# Downloaded files in the folder "omim_searches"


file_list <- list.files(path = "data/omim_searches")

omim_searches = list()

for (i in file_list){
  title = str_replace(i, "_OMIM-Entry-Retrieval.txt", "")
  df = read_delim(paste0("data/omim_searches/", i),
                  "\t", escape_double = FALSE, trim_ws = TRUE, skip = 3) %>%
    mutate(term = title)
  omim_searches[[title]] = df
}

omim_searches = do.call(rbind, omim_searches) %>% 
  mutate(term = case_when(grepl("intellectual_disability", term) ~ "intellectual_disability", 
                          grepl("motor_development", term) ~ "motor_development",
                          grepl("seizures", term) ~ "seizures",
                          grepl("^ClinicalSynopsis", term) ~ "ClinicalSynopsis",
                          grepl("Neurologic_ClinicalSynopsis", term) ~ "Neurologic_ClinicalSynopsis",
                          TRUE ~ term),
         `MIM Number` = str_replace(`MIM Number`, c("#|%"), "")) %>%
  select(`MIM Number`, term)


## Genes annotated with Clinical Synopsis data for their associated diseases 

# For columns referring to Clinical Synopsis annotations, "TRUE" means that there is at least one phenotypic entry associated with the gene in question that fits the criteria 

omim_morbid_split_annot = omim_morbid_split %>%
  left_join(omim_searches, by = c("phen_MIM_number" = "MIM Number")) %>%
  group_by_at(names(.)[which(names(.) != "term")]) %>%
  summarize(term = paste0(term[which(!is.na(term))], collapse = ", ")) %>%
  mutate(ClinicalSynopsis = case_when(grepl("ClinicalSynopsis", term) ~ TRUE), 
         Neurologic_ClinicalSynopsis = case_when(grepl("Neurologic_ClinicalSynopsis", term) ~ TRUE), 
         intellectual_disability = case_when(grepl("intellectual_disability", term) ~ TRUE),
         seizures = case_when(grepl("seizures", term) ~ TRUE),
         motor_development = case_when(grepl("motor_development", term) ~ TRUE),
         language_development = case_when(grepl("language_development", term) ~ TRUE),
         spasticity = case_when(grepl("spasticity", term) ~ TRUE),
         ataxia = case_when(grepl("ataxia", term) ~ TRUE)) %>%
  select(-term)


# The phenotypes with brackets [ ], braces { } and  question mark (?) were moved to a different column ("phen_additional")
# The same was done to their respective MIM numbers and mode of inheritence 

omim_morbid_split_additional = omim_morbid_split_annot %>%
  rowwise() %>%
  mutate(phen_additional = case_when(str_detect(phen, "^\\{") ~ phen,
                                     str_detect(phen, "^\\?") ~ phen, 
                                     str_detect(phen, "^\\[") ~ phen),
         phen_MIM_number_additional = case_when(str_detect(phen, "^\\{") ~ phen_MIM_number,
                                                str_detect(phen, "^\\?") ~ phen_MIM_number, 
                                                str_detect(phen, "^\\[") ~ phen_MIM_number),
         inheritance_additional = case_when(str_detect(phen, "^\\{") ~ inheritance,
                                            str_detect(phen, "^\\?") ~ inheritance, 
                                            str_detect(phen, "^\\[") ~ inheritance)) %>%
  mutate(inheritance = case_when(str_detect(phen, "^\\{") ~ NA_character_,
                                 str_detect(phen, "^\\?") ~ NA_character_,
                                 str_detect(phen, "^\\[") ~ NA_character_,
                                 TRUE ~ inheritance),
         phen_MIM_number = case_when(str_detect(phen, "^\\{") ~ NA_character_,
                                     str_detect(phen, "^\\?") ~ NA_character_,
                                     str_detect(phen, "^\\[") ~ NA_character_,
                                     TRUE ~ phen_MIM_number),
         phen = case_when(str_detect(phen, "^\\{") ~ NA_character_,
                          str_detect(phen, "^\\?") ~ NA_character_,
                          str_detect(phen, "^\\[") ~ NA_character_,
                          TRUE ~ phen))


# Summarize the two groups of phenotypes separately

omim_morbid_split_additional_1 = omim_morbid_split_additional %>%
  select(-c(phen_additional, phen_MIM_number_additional, inheritance_additional)) %>%
  filter(!is.na(phen)) %>%
  mutate(inheritance = case_when(is.na(inheritance) ~ "",
                                 TRUE ~ inheritance)) %>%
  group_by(Chromosome, `Genomic Position Start`, `Genomic Position End`, `Approved Symbol`,
           `Cyto Location`, `Computed Cyto Location`, `MIM Number`, `Gene Symbols`,
           `Gene Name`, `Entrez Gene ID`, `Ensembl Gene ID`, Comments, `Mouse Gene Symbol/ID`) %>%
  summarize(phen = paste0(phen[which(!is.na(phen))], collapse = " // "),
            phen_MIM_number = paste0(phen_MIM_number[which(!is.na(phen_MIM_number))], collapse = " // "),
            inheritance = paste0(inheritance, collapse = " // ")) %>%
  ungroup() %>%
  mutate(phen = case_when(phen == "" ~ NA_character_,
                          TRUE ~ phen),
         phen_MIM_number = case_when(phen_MIM_number == "" ~ NA_character_,
                                     TRUE ~ phen_MIM_number),
         inheritance = case_when(inheritance == "" ~ NA_character_,
                                 TRUE ~ inheritance))

omim_morbid_split_additional_2 = omim_morbid_split_additional %>%
  select(-c(phen, phen_MIM_number, inheritance)) %>%
  filter(!is.na(phen_additional)) %>%
  mutate(inheritance_additional = case_when(is.na(inheritance_additional) ~ "",
                                            TRUE ~ inheritance_additional),
         phen_MIM_number_additional = case_when(is.na(phen_MIM_number_additional) ~ "",
                                                TRUE ~ phen_MIM_number_additional)) %>%
  group_by(Chromosome, `Genomic Position Start`, `Genomic Position End`, `Approved Symbol`,
           `Cyto Location`, `Computed Cyto Location`, `MIM Number`, `Gene Symbols`,
           `Gene Name`, `Entrez Gene ID`, `Ensembl Gene ID`, Comments, `Mouse Gene Symbol/ID`) %>%
  summarize(phen_additional = paste0(phen_additional[which(!is.na(phen_additional))], collapse = " // "),
            phen_MIM_number_additional = paste0(phen_MIM_number_additional, collapse = " // "),
            inheritance_additional = paste0(inheritance_additional, collapse = " // ")) %>%
  ungroup() %>%
  mutate(phen_additional = case_when(phen_additional == "" ~ NA_character_,
                                     TRUE ~ phen_additional),
         phen_MIM_number_additional = case_when(phen_MIM_number_additional == "" ~ NA_character_,
                                                TRUE ~ phen_MIM_number_additional),
         inheritance_additional = case_when(inheritance_additional == "" ~ NA_character_,
                                            TRUE ~ inheritance_additional))


# Join the data keeping the initial order

omim_morbid_modified = omim_morbid_split_additional %>%
  select(-c(phen, phen_MIM_number, inheritance,
            phen_additional, phen_MIM_number_additional, inheritance_additional)) %>%
  group_by(Chromosome, `Genomic Position Start`, `Genomic Position End`, `Approved Symbol`,
           `Cyto Location`, `Computed Cyto Location`, `MIM Number`, `Gene Symbols`,
           `Gene Name`, `Entrez Gene ID`, `Ensembl Gene ID`, Comments, `Mouse Gene Symbol/ID`) %>%
  summarize(Phenotypes = paste0(Phenotypes[which(!is.na(Phenotypes))], collapse = " // "),
            ClinicalSynopsis = any(ClinicalSynopsis),
            Neurologic_ClinicalSynopsis = any(Neurologic_ClinicalSynopsis),
            intellectual_disability = any(intellectual_disability), 
            seizures = any(seizures),
            motor_development = any(motor_development),
            language_development = any(language_development),
            spasticity = any(spasticity),
            ataxia = any(ataxia)) %>%
  ungroup() %>%
  left_join(omim_morbid_split_additional_1) %>%
  left_join(omim_morbid_split_additional_2)


### Merge HGNC and OMIM annotations
############################################################################################################################

df_allChr_allGenes = HGNC_symbols_allChr %>%
  left_join(omim_genes_hg38_annot %>%
              filter(!is.na(`Approved Symbol`)) %>%
              select(`Approved Symbol`, `MIM Number`, morbid) %>%
              group_by(`Approved Symbol`) %>%
              summarize(MIM_Number = paste0(unique(`MIM Number`), collapse = ", "),
                        morbid = any(morbid)) %>%
              ungroup(),
            by = c(`Approved symbol` = "Approved Symbol")) %>%
  left_join(omim_morbid_modified %>%
              select(`Approved Symbol`,
                     phen,
                     phen_MIM_number,
                     inheritance,
                     phen_additional,
                     phen_MIM_number_additional,
                     inheritance_additional,
                     ClinicalSynopsis,
                     Neurologic_ClinicalSynopsis,
                     intellectual_disability,
                     seizures,
                     motor_development,
                     language_development,
                     spasticity,
                     ataxia) %>%
              ungroup(), 
            by = c(`Approved symbol` = "Approved Symbol")) %>%
  unique()


# 18 morbid genes with approved symbol that are non-protein-coding
# Following analyses focused on "protein-coding genes" with "approved" status that are present in the reference assembly  
# Merged HGNC+OMIM data on approved protein-coding genes from all chromsomes was saved in the file "HGNCgenes_OMIMannotated.csv"

df_allChr = df_allChr_allGenes %>%
  filter(all(Status == "Approved", `Locus group` == "protein-coding gene", !grepl("not on reference assembly", Chromosome))) %>%
  unique() %>%
  mutate(morbid = case_when((morbid %in% c(T, F)) ~ morbid,
                            TRUE ~ FALSE),
         gene_category = case_when(!is.na(phen) & ClinicalSynopsis == TRUE ~ "confirmed +CS",
                                   !is.na(phen) & is.na(ClinicalSynopsis) ~ "confirmed -CS",
                                   !is.na(phen_additional) & ClinicalSynopsis == TRUE ~ "PMT +CS",
                                   !is.na(phen_additional) & is.na(ClinicalSynopsis) ~ "PMT -CS",
                                   TRUE ~ "no disorder"))

df_allChr = df_allChr %>%
  ungroup() %>%
  filter(`Approved symbol` != "MED14OS")

write_csv(df_allChr, "results/HGNC_OMIM_all_chr/HGNCgenes_OMIMannotated.csv", na = "")

