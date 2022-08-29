# author: "Elsa Leitao"

library(readr)
library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(matrixStats)
library(tibble)
library(purrr)
library(seqinr)
library(biomaRt)
library(GenomicScores)
library(phastCons100way.UCSC.hg38)


### Download Brainspan and GTEx data  
############################################################################################################################

## Brainspan

# https://www.brainspan.org/static/download.html
# RNA-Seq Gencode v10 summarized to genes
# gene_matrix_csv.zip
# unzip to raw/Brainspan_genes_matrix_csv

## GTEx

# https://www.gtexportal.org/home/datasets)
# GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
# save it in raw/


### Functions
############################################################################################################################
# In some annotation tables, the gene names don't follow the HUGO nomenclature
# Whenever that happens, I convert gene names to HGNC approved symbols 

approved_symbols <- read_delim("raw/HGNC_symbols_allChr.txt",
                               "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  filter(Status == "Approved" & `Locus group` == "protein-coding gene" & !grepl("not on reference assembly", Chromosome))


## Functions to convert gene symbols to HGNC approved symbols

convert_genes_to_approved = function(vector){
  approved = vector[which(vector %in% approved_symbols$`Approved symbol`)]
  not_approved = vector[which(!vector %in% approved_symbols$`Approved symbol`)]
  
  print(length(not_approved))
  if(length(not_approved) != 0){
    counter = 0
    converted = list()
    unconverted = list()
    for(gene in not_approved){
      counter = counter + 1
      print(counter)
      # print(gene)
      rows = approved_symbols %>%
        filter_all(any_vars(grepl(paste0(gene, "$|", gene, ","),., ignore.case = TRUE)))
      if(nrow(rows) == 0){
        unconverted[[gene]] = gene
      } else {
        symbol = rows$`Approved symbol`
        df = data.frame(approved_symbol = symbol,
                        previous_symbol = gene)
        converted[[gene]] = df
      }
    }
    
    converted = do.call(rbind, converted) %>%
      mutate(status = "converted")
    
    approved_df = data.frame(approved_symbol = approved,
                             previous_symbol = approved,
                             status = "approved") %>%
      full_join(converted)
    
    if(length(unconverted) != 0){
      unconverted = do.call(rbind, unconverted) %>%
        data.frame() %>%
        pull() %>%
        as.character()
      return(list(approved_df, unconverted))
      
    } else {
      return(list(approved_df, "all genes with approved symbol"))
    }
  } else {
    return("list of genes contains only approved_symbols")
  }
}

convert_single_gene_to_approved = function(string){
  #print(string)
  if(string %in% approved_symbols$`Approved symbol`){
    return(as.character(string))
  } else {
    symbols = approved_symbols %>% 
      filter_all(any_vars(grepl(paste0(string, "$|", string, ","),., ignore.case = TRUE))) %>%
      pull(`Approved symbol`) %>%
      paste0(collapse = ", ")
    return(symbols)
  }
}


### Import HGNC genes with OMIM annotation from all chromosomes
############################################################################################################################

# file generated in "HGNC_OMIM_genes_all_chr.R"

HGNC_OMIM <- read_csv("results/HGNC_OMIM_all_chr/HGNCgenes_OMIMannotated.csv") %>%
  filter(Status == "Approved") %>%
  mutate(monogenic_disorder = case_when(gene_category %in% c("confirmed +CS", "confirmed -CS") ~ "confirmed",
                                        gene_category %in% c("PMT +CS", "PMT -CS") ~ "PMT",
                                        gene_category == "no disorder" ~ "no disorder")) %>%
  ungroup()


### gnomAD annotations 
############################################################################################################################

# file from "https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
# In OMIM tables, the gene names are HGNC approved symbols, but in gnomAD they aren't
# Before joining data, I converted gnomAD gene names to HGNC approved symbols. However, some of them had no approved symbols.

# I kept the raw data columns 

# In the gnomad table there where a few genes occuring more than once (multiple transcripts)
# I kept the occurrence with smaller oef_lof_upper (= LOEUF).


gnomad_preslice_pre = read.delim("raw/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz") 


## Retrieve Ensembl gene ID based on Ensembl transcript ID
# biomaRt (use hg19 for gnomAD v2.1.1)

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://grch37.ensembl.org")


# Bring coding exons coordinates from canonical transcripts  

geneIDs = getBM(attributes = c("ensembl_transcript_id",
                               "ensembl_gene_id"),
                filters = c("ensembl_transcript_id"),
                values = list(gnomad_preslice_pre$transcript),
                mart = ensembl) 


gnomad_preslice = gnomad_preslice_pre %>%
  left_join(geneIDs, by = c(transcript = "ensembl_transcript_id")) %>%
  dplyr::select(ensembl_gene_id, gene:exp_hom_lof, -constraint_flag, -oe_lof_upper_bin, -oe_lof_upper_bin_6, num_coding_exons, gene_length) %>%
  unique()


### GTEx annotations 
############################################################################################################################

## Import GTEX data

# GTEx the file "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz" 
# containing the median gene-level TPM by tissue (54 tissues) from all chromosomes

# (long calculation time (start))

GTEx <- read_table2("raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",
    skip = 2) 

# (long calculation time (end))


## Restrict to protein coding genes in the HGNC_OMIM table based on the Ensembl ID

# 4 genes don't have the same Ensembl ID in the 2 columns
# use the ensemblID that is present in the GTEx table

table(HGNC_OMIM$`Ensembl ID(supplied by Ensembl)` == HGNC_OMIM$`Ensembl gene ID`)
HGNC_OMIM %>%
  filter(HGNC_OMIM$`Ensembl ID(supplied by Ensembl)` != HGNC_OMIM$`Ensembl gene ID`) %>%
  select(`Approved symbol`, `Ensembl gene ID`, `Ensembl ID(supplied by Ensembl)`, `Previous symbols`, `Alias symbols`)

HGNC_OMIM = HGNC_OMIM %>%
  mutate(ensemblID = case_when(`Ensembl ID(supplied by Ensembl)` %in% c("ENSG00000276085", "ENSG00000205045") ~ `Ensembl ID(supplied by Ensembl)`,
                               `Ensembl gene ID` == "ENSG00000286065" ~ `Ensembl gene ID`,
                               TRUE ~ `Ensembl gene ID`))

genes = data.frame(Name = GTEx$Name) %>%
  rowwise() %>%
  mutate(Name2 = str_split(Name, "[.]") %>%
           unlist() %>%
           .[1])

GTExgenes_in_HGNCOMIM = genes %>%
  filter(Name2 %in% HGNC_OMIM$ensemblID)

GTEx_filt = GTEx %>%
  filter(Name %in% GTExgenes_in_HGNCOMIM$Name)


## Import GTEx sample annotations

GTEx_samples_annot <- read_delim("raw/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  select(SUBJID, SEX)

GTEx_sample <- read_delim("raw/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  filter(SMAFRZE == "RNASEQ") %>%
  select(SAMPID, SMTSD) %>%
  rowwise() %>%
  mutate(SUBJID = str_split(SAMPID, "-") %>% unlist() %>% .[1:2] %>% paste0(collapse = "-")) %>%
  left_join(GTEx_samples_annot)


## Group tissues by similar types 

GTEx_sample_grouped = GTEx_sample %>%
  mutate(tissue = case_when(SMTSD == "Nerve - Tibial" ~ "Nerve",
                            SMTSD %in% c("Brain - Cerebellar Hemisphere", "Brain - Cerebellum") ~ "Cerebellar_tissues",
                            SMTSD %in% c("Brain - Amygdala",
                                         "Brain - Anterior cingulate cortex (BA24)",
                                         "Brain - Caudate (basal ganglia)",
                                         "Brain - Cortex",
                                         "Brain - Frontal Cortex (BA9)",
                                         "Brain - Hippocampus",
                                         "Brain - Hypothalamus",
                                         "Brain - Nucleus accumbens (basal ganglia)",
                                         "Brain - Putamen (basal ganglia)",
                                         "Brain - Spinal cord (cervical c-1)",
                                         "Brain - Substantia nigra") ~ "Other_Brain_tissues",
                            TRUE ~ "Other_tissues"))


## Expression separating females and males

expression_groups = list()

for (i in 1:nrow(GTEx_groups)){
  Tissue = GTEx_groups$tissue[i]
  Sex = GTEx_groups$SEX[i]
  colNameMean = paste0(Tissue, "_mean_", as.character(Sex))
  colNameVar = paste0(Tissue, "_var_", as.character(Sex))
  
  samples = GTEx_sample_grouped %>%
    filter(tissue == Tissue & SEX == Sex) %>%
    pull(SAMPID)
  
  data = GTEx_filt %>%
    select(Name, Description, one_of(samples)) %>%
    mutate(mean = rowMeans(.[, 3:ncol(.)], na.rm = TRUE),
           var = rowVars(as.matrix(.[, 3:ncol(.)]), na.rm = TRUE)) %>%
    select(Name, Description, mean, var) %>%
    rename(!!colNameMean := mean,
           !!colNameVar := var)
  
  expression_groups[[i]] = data
}

expression_groups2 = purrr::reduce(expression_groups, full_join, by = c("Name", "Description"))


### tau annotations (using GTEx data)
############################################################################################################################

# "tau" as a measure of how tissue specific the expression of a gene is:
# - tau < 0.6 (broadly expressed)  
# - tau > 0.6 (expressed in a restricted number of tissues)  


## Calculate median of median expression in cerebellar tissues (2 regions) and in the other brain tissues (11 regions)
# Discard the expression in individual brain tissues
# Group tissues by similar types only for cerebellum and other brain regions (43 tissue groups)

GTEx_sample_grouped_tau = GTEx_sample %>%
  mutate(tissue = case_when(SMTSD %in% c("Brain - Cerebellar Hemisphere", "Brain - Cerebellum") ~ "Cerebellar_tissues",
                            SMTSD %in% c("Brain - Amygdala",
                                         "Brain - Anterior cingulate cortex (BA24)",
                                         "Brain - Caudate (basal ganglia)",
                                         "Brain - Cortex",
                                         "Brain - Frontal Cortex (BA9)",
                                         "Brain - Hippocampus",
                                         "Brain - Hypothalamus",
                                         "Brain - Nucleus accumbens (basal ganglia)",
                                         "Brain - Putamen (basal ganglia)",
                                         "Brain - Spinal cord (cervical c-1)",
                                         "Brain - Substantia nigra") ~ "Other_Brain_tissues",
                            TRUE ~ SMTSD))


## Calculate tau (all individuals together)

# Median TPM Expression

expression_groups_tau = list()

for (i in 1:nrow(GTEx_groups_tau)){
  Tissue = GTEx_groups_tau$tissue[i]
  colNameMedian = paste0(Tissue, "_median")
  
  samples = GTEx_sample_grouped_tau %>%
    filter(tissue == Tissue) %>%
    pull(SAMPID)
  
  data = GTEx_filt %>%
    select(Name, Description, one_of(samples)) %>%
    mutate(median = rowMedians(as.matrix(.[, 3:ncol(.)]), na.rm = TRUE)) %>%
    select(Name, Description, median) %>%
    rename(!!colNameMedian := median)
  
  expression_groups_tau[[i]] = data
}

expression_groups2_tau = purrr::reduce(expression_groups_tau, full_join, by = c("Name", "Description"))


# Calculate log2(TPM+1)

expression_groups2_tau_log <- expression_groups2_tau %>%
  mutate_if(is.double, funs(log2(.+1)))


# Calculate the maximal expression value

expression_groups2_tau_log = expression_groups2_tau_log %>%
  gather(tissue, value, 3:ncol(.)) %>%
  group_by(Name) %>%
  summarize(max_expression = max(value)) %>%
  right_join(expression_groups2_tau_log, by = "Name") %>%
  dplyr::select(-max_expression, max_expression) 


# Calculate expression profile component normalized by the maximal component value

expression_groups2_tau_log_normlog = expression_groups2_tau_log %>%
  mutate_if(is.double, funs(./max_expression)) %>%
  dplyr::select(-max_expression)


# Calculate tau

expression_groups2_tau_log_normlog_1 = expression_groups2_tau_log_normlog %>%
  mutate_if(is.double, funs(1-.)) 

expression_groups2_tau_log_normlog_1_tau = expression_groups2_tau_log_normlog_1 %>%
  gather(tissue, value, 3:ncol(.)) %>%
  group_by(Name) %>%
  summarize(tau = sum(value) / n()) %>%
  right_join(expression_groups2_tau_log_normlog_1, by = "Name") %>%
  dplyr::select(Name, Description, tau)



## Calculate tau (female and male separate)

# Median TPM Expression

expression_groups_tau_sex = list()

for (i in 1:nrow(GTEx_groups_tau_sex)){
  Tissue = GTEx_groups_tau_sex$tissue[i]
  Sex = GTEx_groups_tau_sex$SEX[i]
  colNameMedian = paste0(Tissue, "_median_", as.character(Sex))
  
  samples = GTEx_sample_grouped_tau %>%
    filter(tissue == Tissue & SEX == Sex) %>%
    pull(SAMPID)
  
  data = GTEx_filt %>%
    select(Name, Description, one_of(samples)) %>%
    mutate(median = rowMedians(as.matrix(.[, 3:ncol(.)]), na.rm = TRUE)) %>%
    select(Name, Description, median) %>%
    rename(!!colNameMedian := median)
  
  expression_groups_tau_sex[[i]] = data
}

expression_groups2_tau_sex = purrr::reduce(expression_groups_tau_sex, full_join, by = c("Name", "Description"))


# Calculate log2(TPM+1)

expression_groups2_tau_sex_log <- expression_groups2_tau_sex %>%
  mutate_if(is.double, funs(log2(.+1)))


# Calculate the maximal expression value

#- MALES
expression_groups2_tau_Males_log = expression_groups2_tau_sex_log %>%
  select(Name, Description, ends_with("_1")) %>%
  gather(tissue, value, 3:ncol(.)) %>%
  group_by(Name) %>%
  summarize(max_expression = max(value)) %>%
  right_join(expression_groups2_tau_sex_log %>%
               select(Name, Description, ends_with("_1")), by = "Name") %>%
  dplyr::select(-max_expression, max_expression) 

#- FEMALES
expression_groups2_tau_Females_log = expression_groups2_tau_sex_log %>%
  select(Name, Description, ends_with("_2")) %>%
  gather(tissue, value, 3:ncol(.)) %>%
  group_by(Name) %>%
  summarize(max_expression = max(value)) %>%
  right_join(expression_groups2_tau_sex_log %>%
               select(Name, Description, ends_with("_2")), by = "Name") %>%
  dplyr::select(-max_expression, max_expression) 


# Calculate expression profile component normalized by the maximal component value

#- MALES
expression_groups2_tau_Males_log_normlog = expression_groups2_tau_Males_log %>%
  mutate_if(is.double, funs(./max_expression)) %>%
  dplyr::select(-max_expression)

#- FEMALES
expression_groups2_tau_Females_log_normlog = expression_groups2_tau_Females_log %>%
  mutate_if(is.double, funs(./max_expression)) %>%
  dplyr::select(-max_expression)


# Calculate tau

#- MALES
expression_groups2_tau_Males_log_normlog_1 = expression_groups2_tau_Males_log_normlog %>%
  mutate_if(is.double, funs(1-.)) 

expression_groups2_tau_Males_log_normlog_1_tau = expression_groups2_tau_Males_log_normlog_1 %>%
  gather(tissue, value, 3:ncol(.)) %>%
  group_by(Name) %>%
  summarize(tau_1 = sum(value) / n()) %>%
  right_join(expression_groups2_tau_Males_log_normlog_1, by = "Name") %>%
  dplyr::select(Name, Description, tau_1)


#- FEMALES
expression_groups2_tau_Females_log_normlog_1 = expression_groups2_tau_Females_log_normlog %>%
  mutate_if(is.double, funs(1-.)) 

expression_groups2_tau_Females_log_normlog_1_tau = expression_groups2_tau_Females_log_normlog_1 %>%
  gather(tissue, value, 3:ncol(.)) %>%
  group_by(Name) %>%
  summarize(tau_2 = sum(value) / n()) %>%
  right_join(expression_groups2_tau_Females_log_normlog_1, by = "Name") %>%
  dplyr::select(Name, Description, tau_2)


## tau data together

tau_data = expression_groups2_tau_log_normlog_1_tau %>%
  left_join(expression_groups2_tau_Males_log_normlog_1_tau) %>%
  left_join(expression_groups2_tau_Females_log_normlog_1_tau)


### Import X-inactivation annotations
############################################################################################################################

# file generated in "annotations_chrX_genes.R"

XCI_curated_2 <- read_excel("results/annotations_chrX/curated/XCI_data_curated.xlsx") %>%
  dplyr::select(Gene, XCI) %>%
  unique() %>%
  mutate(XCI = case_when(XCI == "na" ~ NA_character_,
                         TRUE ~ XCI))

                         
### Brainspan annotations
############################################################################################################################

# downloaded from Brainspan the "Developmental Transcriptome Dataset"
# "RNA-Seq Gencode v10 summarized to genes" (https://www.brainspan.org/static/download.html)
# containing RPKM normalized expression values and meta-data. 

# Brainspan data was generated across 13 developmental stages in 8-16 brain structures from a total of 42 individuals of both sexes
# (note: not the same number of individuals nor the same sex distribution across the different developmental stages)
# Additional information available at "https://help.brain-map.org/download/attachments/3506181/Transcriptome_Profiling.pdf?version=1&modificationDate=1382036562736&api=v2"


## Retrieve BRAINSPAN data 

# https://www.brainspan.org/static/download.html
# RNA-Seq Gencode v10 summarized to genes
# gene_matrix_csv.zip
# unzip to Brainspan_genes_matrix_csv


## Import rows_metadata and convert gene names to approved symbols

# (long calculation time (start))

rows_metadata = read_csv("raw/Brainspan_genes_matrix_csv/rows_metadata.csv") %>% 
  rowwise() %>%
  mutate(Gene = convert_single_gene_to_approved(gene_symbol)) %>%
  mutate(Gene = case_when(Gene == "" ~ as.character(gene_symbol),
                          TRUE ~ Gene)) %>%
  select(-gene_symbol)

# (long calculation time (end))


## Correct a few convertions that were missed

for(i in setdiff(HGNC_OMIM$`Approved symbol`, rows_metadata$Gene)){
  rows_metadata = rows_metadata %>%
    mutate(Gene = case_when((grepl(paste0("^", i, "|", ", ", i), Gene) & ensembl_gene_id != "ENSG00000164458") ~ i,
                            TRUE ~ Gene))
}


## Import column data

columns_metadata <- read_csv("raw/Brainspan_genes_matrix_csv/columns_metadata.csv")


## Import expression data in RPKM

expression_matrix <- read_csv("raw/Brainspan_genes_matrix_csv/expression_matrix.csv", 
                              col_names = FALSE)


## Convert RPKM to TPM

RPKMsum <- expression_matrix %>%
  dplyr::select(-X1) %>%
  summarize_all(sum) %>%
  gather(sample, RPKMsum, X2:X525)

expression_matrix = expression_matrix %>%
  gather(sample, RPKM, X2:X525) %>%
  left_join(RPKMsum) %>%
  mutate(TPM = RPKM * 10^6 / RPKMsum)


## Keep data related to genes with HGNC appoved symbols (Row indexes referring to these genes) 

rows_filt = which(rows_metadata$Gene %in% HGNC_OMIM$`Approved symbol`) 

rows_metadata_filt = rows_metadata %>%
  filter(row_num %in% rows_filt)

rows_metadata_curated = rows_metadata %>%
  filter(row_num %in% rows_filt) %>%
  filter(Gene %in% HGNC_OMIM$`Approved symbol` | entrez_id %in% HGNC_OMIM$`NCBI Gene ID(supplied by NCBI)`) %>%
  group_by(Gene) %>%
  arrange(desc(ensembl_gene_id %in% HGNC_OMIM$ensemblID)) %>%
  dplyr::slice(1)

genes_missing = rows_metadata %>%
  filter(rows_metadata$Gene %in% setdiff(rows_metadata_filt$Gene %>% unique, rows_metadata_curated$Gene %>% unique()))

rows_metadata_curated2 = rows_metadata_curated %>%
  full_join(genes_missing)


## Expression data referring to genes in the HGNC_OMIM list

expression_matrix_filt = expression_matrix %>%
  filter(X1 %in% rows_metadata_curated2$row_num) %>%
  left_join(rows_metadata_filt %>% 
              dplyr::select(Gene, row_num), by = c(X1 = "row_num"))


## Check columns_metadata to see which columns from the expression_matrix to summarize

columns_metadata %>%
  mutate(age_group = case_when(str_detect(age, "pcw") ~ "prenatal",
                               age %in% c("4 mos", "10 mos", "1 yrs",  "2 yrs", "3 yrs", "4 yrs") ~ "postnatal_1",
                               age %in% c("8 yrs", "11 yrs", "13 yrs", "15 yrs", "18 yrs", "19 yrs", "21 yrs", "23 yrs", "30 yrs",
                                          "36 yrs", "37 yrs", "40 yrs") ~ "postnatal_2")) %>%
  group_by(age_group, gender) %>%
  summarize(count = n())


# MALE

#- Prenatal: age in pcw
prenatal_cols_M = paste0("X", (which(grepl("pcw", columns_metadata$age) & grepl("M", columns_metadata$gender)) + 1)) 
# (+ 1 because expression_matrix has a first column with the row_num)

#- Post-natal 1: after birth until 4 year-old (inclusive)
postnatal_1_cols_M = paste0("X", which(columns_metadata$age %in% c("4 mos", "10 mos", "1 yrs",  "2 yrs", "3 yrs", "4 yrs") & grepl("M", columns_metadata$gender)) + 1)
# (+ 1 because expression_matrix has a first column with the row_num)

#- Post-natal 2: older or equal to 8 years-old
postnatal_2_cols_M = paste0("X", which(columns_metadata$age %in% c("8 yrs", "11 yrs", "13 yrs", "15 yrs", "18 yrs", "19 yrs", "21 yrs", "23 yrs", "30 yrs",
                                                                   "36 yrs", "37 yrs", "40 yrs") & grepl("M", columns_metadata$gender)) + 1)
# (+ 1 because expression_matrix has a first column with the row_num)


# FEMALE

#- Prenatal: age in pcw
prenatal_cols_F = paste0("X", (which(grepl("pcw", columns_metadata$age) & grepl("F", columns_metadata$gender)) + 1)) 
# (+ 1 because expression_matrix has a first column with the row_num)

#- Post-natal 1: after birth until 4 year-old (inclusive)
postnatal_1_cols_F = paste0("X", which(columns_metadata$age %in% c("4 mos", "10 mos", "1 yrs",  "2 yrs", "3 yrs", "4 yrs") & grepl("F", columns_metadata$gender)) + 1)
# (+ 1 because expression_matrix has a first column with the row_num)

#- Post-natal 2: older or equal to 8 years-old
postnatal_2_cols_F = paste0("X", which(columns_metadata$age %in% c("8 yrs", "11 yrs", "13 yrs", "15 yrs", "18 yrs", "19 yrs", "21 yrs", "23 yrs", "30 yrs",
                                                                   "36 yrs", "37 yrs", "40 yrs") & grepl("F", columns_metadata$gender)) + 1)
# (+ 1 because expression_matrix has a first column with the row_num)


## Simplify expression matrix by averaging on age groups and sex while ignoring and tissue 

expression_by_age = expression_matrix_filt %>%
  dplyr::select(Gene, sample, TPM) %>%
  mutate(age_group = case_when(sample %in% prenatal_cols_M ~ "prenatal_M",
                               sample %in% postnatal_1_cols_M ~ "postnatal_1_M",
                               sample %in% postnatal_2_cols_M ~ "postnatal_2_M",
                               sample %in% prenatal_cols_F ~ "prenatal_F",
                               sample %in% postnatal_1_cols_F ~ "postnatal_1_F",
                               sample %in% postnatal_2_cols_F ~ "postnatal_2_F")) %>%
  mutate(age_group2 = case_when(sample %in% prenatal_cols_M ~ "prenatal_M",
                                sample %in% c(postnatal_1_cols_M, postnatal_2_cols_M) ~ "postnatal_M",
                                sample %in% prenatal_cols_F ~ "prenatal_F",
                                sample %in% c(postnatal_1_cols_F, postnatal_2_cols_F) ~ "postnatal_F"))


## Mean and variance 

brainspan_mean1 = expression_by_age %>%
  group_by(Gene, age_group) %>%
  summarize(mean = mean(TPM),
            variance = var(TPM))

brainspan_mean2 = expression_by_age %>%
  group_by(Gene, age_group2) %>%
  summarize(mean = mean(TPM),
            variance = var(TPM)) %>%
  dplyr::rename(age_group = age_group2) %>%
  filter(age_group %in% c("postnatal_M", "postnatal_F"))

brainspan_mean = brainspan_mean1 %>% 
  rbind(brainspan_mean2) %>%
  arrange(Gene)


## Combine data

brainspan_data_mean <- brainspan_mean %>%
  mutate(column1 = case_when(age_group == "prenatal_M" ~ "mean_Pre_M",
                             age_group == "postnatal_M" ~ "mean_Post_M",
                             age_group == "postnatal_1_M" ~ "mean_Post1_M",
                             age_group == "postnatal_2_M" ~ "mean_Post2_M",
                             age_group == "prenatal_F" ~ "mean_Pre_F",
                             age_group == "postnatal_F" ~ "mean_Post_F",
                             age_group == "postnatal_1_F" ~ "mean_Post1_F",
                             age_group == "postnatal_2_F" ~ "mean_Post2_F")) %>%
  dplyr::select(-age_group, -variance) %>%
  spread(column1, mean) 

brainspan_data_var <- brainspan_mean %>%
  mutate(column1 = case_when(age_group == "prenatal_M" ~ "var_Pre_M",
                             age_group == "postnatal_M" ~ "var_Post_M",
                             age_group == "postnatal_1_M" ~ "var_Post1_M",
                             age_group == "postnatal_2_M" ~ "var_Post2_M",
                             age_group == "prenatal_F" ~ "var_Pre_F",
                             age_group == "postnatal_F" ~ "var_Post_F",
                             age_group == "postnatal_1_F" ~ "var_Post1_F",
                             age_group == "postnatal_2_F" ~ "var_Post2_F")) %>%
  dplyr::select(-age_group, -mean) %>%
  spread(column1, variance) 

brainspan_data = full_join(brainspan_data_mean, brainspan_data_var)


### Promoter CpG density
############################################################################################################################

# Transcripts retrieved from Ensembl via BioMart using the filters/attributes below :

# Ensembl Genes 103  
# Dataset:  
# - Human genes (GRCh38.p13)  
# Filters:  
# - Transcript type: protein_coding  
# Attributes:  
# - Gene stable ID  
# - Transcript stable ID  
# - Transcript stable ID version  
# - HGNC ID  
# - HGNC symbol  
# - Chromosome/scaffold name  
# - Gene start (bp)  
# - Gene end (bp)  
# - Transcript start (bp)  
# - Transcript end (bp)  
# - Transcription start site (TSS)  
# - GENCODE basic annotation  
# - APPRIS annotation  
# - RefSeq match transcript (MANE Select)  
# saved as "mart_export.txt"

# For restricting to canonical transcripts I considered MANE and APPRIS annotations:  
# - MANE Select transcripts (Matched Annotation between NCBI and EBI) were independently identified by both Ensembl and NCBI as the most biologically relevant.   
# - APPRIS is a system to annotate alternatively spliced transcripts based on a range of computational methods.  
# More information about the transcript quality tags can be found in https://www.ensembl.org/info/genome/genebuild/transcript_quality_tags.html.

# Note: for one case I had to use the "Ensembl ID(supplied by Ensembl)" and not the HGNC curated version because it was wrong.  
# The ratio observed/expected (Obs/Exp or O/E) CpG was calculated as follows: Obs/Exp CpG = (Number of CpG / (Number of C x Number of G)) x N, where N is the total number of nucleotides in the sequence being analysed.  

ensembl_transcripts <- read_delim("raw/mart_export.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE,
                                  col_types = cols(`RefSeq match transcript (MANE Select)` = col_character())) %>%
  mutate(length = `Transcript end (bp)` - `Transcript start (bp)`) %>%
  filter(`Chromosome/scaffold name` %in% c(as.character(seq(1:22)), "X", "Y"))


## MANE and APPRIS annotations
# - MANE Select transcripts (Matched Annotation between NCBI and EBI) were independently identified by both Ensembl and NCBI as the most biologically relevant.   
# - APPRIS is a system to annotate alternatively spliced transcripts based on a range of computational methods.
MANE_transcripts = ensembl_transcripts %>% 
  filter(!is.na(`RefSeq match transcript (MANE Select)`) & `Gene stable ID` %in% HGNC_OMIM$ensemblID)

APPRIS_transcripts = ensembl_transcripts %>% 
  filter(!`Gene stable ID` %in% MANE_transcripts$`Gene stable ID` & !is.na(`APPRIS annotation`) & `Gene stable ID` %in% HGNC_OMIM$ensemblID) %>%
  group_by(`Gene stable ID`) %>%
  arrange(desc(length)) %>%
  dplyr::slice(1)

other_transcripts = ensembl_transcripts %>%
  filter(`Gene stable ID` %in% HGNC_OMIM$ensemblID &
           !`Gene stable ID` %in% MANE_transcripts$`Gene stable ID` &
           !`Gene stable ID` %in% APPRIS_transcripts$`Gene stable ID`) %>%
  group_by(`Gene stable ID`) %>%
  arrange(desc(length)) %>%
  dplyr::slice(1)


## Canonical transcripts

# The rules to consider canonical transcripts followed this order:  
# 1 - MANE transcript if it exists (vast majority. approx. 80% genes)  
# 2 - longest transcript with APPRIS annotation (approx. 20% genes)  
# 3 - existing transcript (just in one case)  

canonical_transcripts = MANE_transcripts %>%
  full_join(APPRIS_transcripts) %>%
  full_join(other_transcripts) 


## Coordinates promoters canonical transcripts
# 4kb regions surrounding the TSS (+-2kb)

canonical_TSS = canonical_transcripts %>%
  dplyr::select(`Chromosome/scaffold name`, `Transcription start site (TSS)`) %>%
  mutate(start = `Transcription start site (TSS)`- 2000,
         end = `Transcription start site (TSS)`+ 2000) %>%
  dplyr::select(`Chromosome/scaffold name`, start, end) %>%
  write_tsv( "results/annotations_all_chr_ML/promoters.bed", col_names = FALSE)


## Get promoter sequences
# "results/annotations_all_chr_ML/promoters.bed" uploaded as a "Custom Track" in UCSC and obtained the sequences through "Table Browser"
# saved as "raw/promoter_seqs.txt"

promoter_seqs = read.fasta("raw/promoter_seqs.txt", as.string = TRUE)


## CpG density

CpG_density_list = list()
for (i in 1:length(promoter_seqs)){
  seq = promoter_seqs[[i]][1]
  coordinates = attributes(promoter_seqs[[i]])$Annot %>% 
    str_split(" ") %>% 
    unlist() %>% 
    .[2] %>% 
    str_replace("range=", "")
  length = str_length(seq)
  C =  seq %>%
    str_count("c")
  G = seq %>%
    str_count("g")
  CG = seq %>%
    str_count("cg")
  density = (CG / (C*G))*length
  CpG_density_list[[i]] = data.frame(promoter = coordinates,
                                     CpG_density = density)
}

CpG_density_list = do.call(rbind, CpG_density_list) %>%
  unique()

canonical_transcripts_2 = canonical_transcripts %>%
  mutate(tmp1 = paste0("chr", `Chromosome/scaffold name`),
         start = `Transcription start site (TSS)`- 1999,
         end = `Transcription start site (TSS)`+ 2000) %>%
  unite(tmp2, tmp1:start, sep = ":") %>%
  unite(promoter, tmp2:end, sep = "-") %>%
  left_join(CpG_density_list) %>%
  dplyr::select(`Gene stable ID`, CpG_density)


### Promoter and exon-conservation across species
############################################################################################################################

## Use conservation scores from phastCons100way.UCSC.hg38
# PhastCons score ranges from 0 to 1 and represents the probability that a given nucleotide is conserved

phast <- phastCons100way.UCSC.hg38
class(phast)


## Exons

# biomaRt 
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")


# Bring coding exons coordinates from canonical transcripts  
exons_coordinates = getBM(attributes = c("ensembl_transcript_id",
                                         "chromosome_name",
                                         "exon_chrom_start",
                                         "exon_chrom_end",
                                         "cds_start",
                                         "cds_end",
                                         "genomic_coding_start",
                                         "genomic_coding_end"),
                          filters = c("ensembl_transcript_id"),
                          values = list(canonical_transcripts$`Transcript stable ID`),
                          mart = ensembl) %>%
  mutate(chromosome_name = paste0("chr", chromosome_name),
         all_coding = case_when((cds_end - cds_start + 1) == (exon_chrom_end - exon_chrom_start + 1) ~ TRUE,
                                TRUE ~ FALSE)) %>%
  na.omit() %>% # removes non-coding exons
  group_by(ensembl_transcript_id) %>%
  mutate(start = genomic_coding_start,
         end = genomic_coding_end) %>%
  left_join(canonical_transcripts %>% 
              dplyr::select(`HGNC symbol`, `Gene stable ID`, `Transcript stable ID`), 
            by = c(ensembl_transcript_id = "Transcript stable ID"))


## Scores from coding exons

# (long calculation time (start))

exon_scores = list()
for(i in unique(exons_coordinates$ensembl_transcript_id)){

  df_exons = exons_coordinates %>%
    filter(ensembl_transcript_id == i)

  name = df_exons$`HGNC symbol` %>% unique()

  nucleotide_scores = list()
  for(j in 1:nrow(df_exons)){
    chr = df_exons$chromosome_name[j]
    s = df_exons$start[j]
    e = df_exons$end[j]
    scores = gscores(phast, GRanges(seqnames=chr, IRanges(start=s:e, width=1))) %>%
      data.frame()
    nucleotide_scores[[j]] = scores
  }
  nucleotide_scores = do.call(rbind, nucleotide_scores)
  average = mean(nucleotide_scores$default, na.rm = TRUE)

  data = data.frame(transcript_ID = i,
                    HGNC_symbol = name,
                    exon_score = average)
  exon_scores[[i]] = data
}

exon_scores_2 = do.call(rbind, exon_scores)


## Promoters

promoter_coordinates = canonical_transcripts %>% # from the CpG density calculations 
  dplyr::select(`Transcript stable ID`, `HGNC symbol`, `Transcription start site (TSS)`, `Chromosome/scaffold name`) %>%
  mutate(start = `Transcription start site (TSS)`- 2000,
         end = `Transcription start site (TSS)`+ 2000)

promoter_scores = list()
for(i in 1:nrow(promoter_coordinates)){
  chr = promoter_coordinates$`Chromosome/scaffold name`[i]
  s = promoter_coordinates$start[i]
  e = promoter_coordinates$end[i]
  scores = gscores(phast, GRanges(seqnames=chr, IRanges(start=s:e, width=1))) %>%
      data.frame()
  average = mean(scores$default, na.rm=TRUE)

  data = data.frame(transcript_ID = promoter_coordinates$`Transcript stable ID`[i],
                    HGNC_symbol = promoter_coordinates$`HGNC symbol`[i],
                    promoter_score = average)
  promoter_scores[[i]] = data
}

promoter_scores_2 = do.call(rbind, promoter_scores)

# (long calculation time (end))


## Both conservation scores

conservation_scores = exon_scores_2 %>%
  left_join(promoter_scores_2)


### CDS_length
############################################################################################################################

length = exons_coordinates %>%
  mutate(CDS_length = end - start + 1) %>%
  group_by(ensembl_transcript_id) %>%
  summarize(CDS_length = sum(CDS_length))


## Data together

conservation_scores_CDS_length = conservation_scores %>%
  left_join(length, by = c(transcript_ID = "ensembl_transcript_id"))


### Paralogues from Ensembl
############################################################################################################################

# "Paralogues are defined in Ensembl as genes for which the most common ancestor node is a duplication event
# These ancestral duplications are represented by red nodes in the gene trees.
# The table shows the taxonomic level of the ancestor duplication node,
# the Ensembl gene ID and name, the location of the paralogue, 
# and the percent of identical amino acids in the paralogue compared with the gene of interest (Target %ID)
# The identity of the gene of interest when compared with the paralogue is the query %ID."


## Bring paralogues genes using biomaRt  

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")

paralogues1 = getBM(attributes = c("ensembl_gene_id",
                                   "hsapiens_paralog_associated_gene_name",
                                   "hsapiens_paralog_chromosome",
                                   "hsapiens_paralog_orthology_type",
                                   "hsapiens_paralog_subtype",
                                   "hsapiens_paralog_perc_id",
                                   "hsapiens_paralog_perc_id_r1"),
                    filters = c("ensembl_gene_id"),
                    values = list(canonical_transcripts$`Gene stable ID`[1:5000]),
                    mart = ensembl) %>%
  left_join(canonical_transcripts %>%
              dplyr::select(`HGNC symbol`, `Gene stable ID`), by = c(ensembl_gene_id = "Gene stable ID"))

paralogues2 = getBM(attributes = c("ensembl_gene_id",
                                   "hsapiens_paralog_associated_gene_name",
                                   "hsapiens_paralog_chromosome",
                                   "hsapiens_paralog_orthology_type",
                                   "hsapiens_paralog_subtype",
                                   "hsapiens_paralog_perc_id",
                                   "hsapiens_paralog_perc_id_r1"),
                    filters = c("ensembl_gene_id"),
                    values = list(canonical_transcripts$`Gene stable ID`[5001:10000]),
                    mart = ensembl) %>%
  left_join(canonical_transcripts %>%
              dplyr::select(`HGNC symbol`, `Gene stable ID`), by = c(ensembl_gene_id = "Gene stable ID"))

paralogues3 = getBM(attributes = c("ensembl_gene_id",
                                   "hsapiens_paralog_associated_gene_name",
                                   "hsapiens_paralog_chromosome",
                                   "hsapiens_paralog_orthology_type",
                                   "hsapiens_paralog_subtype",
                                   "hsapiens_paralog_perc_id",
                                   "hsapiens_paralog_perc_id_r1"),
                    filters = c("ensembl_gene_id"),
                    values = list(canonical_transcripts$`Gene stable ID`[10001:12000]),
                    mart = ensembl) %>%
  left_join(canonical_transcripts %>%
              dplyr::select(`HGNC symbol`, `Gene stable ID`), by = c(ensembl_gene_id = "Gene stable ID"))

paralogues4 = getBM(attributes = c("ensembl_gene_id",
                                   "hsapiens_paralog_associated_gene_name",
                                   "hsapiens_paralog_chromosome",
                                   "hsapiens_paralog_orthology_type",
                                   "hsapiens_paralog_subtype",
                                   "hsapiens_paralog_perc_id",
                                   "hsapiens_paralog_perc_id_r1"),
                    filters = c("ensembl_gene_id"),
                    values = list(canonical_transcripts$`Gene stable ID`[12001:15000]),
                    mart = ensembl) %>%
  left_join(canonical_transcripts %>%
              dplyr::select(`HGNC symbol`, `Gene stable ID`), by = c(ensembl_gene_id = "Gene stable ID"))

paralogues5 = getBM(attributes = c("ensembl_gene_id",
                                   "hsapiens_paralog_associated_gene_name",
                                   "hsapiens_paralog_chromosome",
                                   "hsapiens_paralog_orthology_type",
                                   "hsapiens_paralog_subtype",
                                   "hsapiens_paralog_perc_id",
                                   "hsapiens_paralog_perc_id_r1"),
                    filters = c("ensembl_gene_id"),
                    values = list(canonical_transcripts$`Gene stable ID`[15001:length(canonical_transcripts$`Gene stable ID`)]),
                    mart = ensembl) %>%
  left_join(canonical_transcripts %>%
              dplyr::select(`HGNC symbol`, `Gene stable ID`), by = c(ensembl_gene_id = "Gene stable ID"))

paralogues = rbind(paralogues1, paralogues2) %>%
  rbind(paralogues3) %>%
  rbind(paralogues4) %>%
  rbind(paralogues5)


## Calculate percentiles of the Target% and Query%

Target_98_perc = quantile(paralogues$hsapiens_paralog_perc_id, probs = seq(0, 1, 0.01), na.rm = TRUE)[99] # 98% percentile Target%

Query_98_perc = quantile(paralogues$hsapiens_paralog_perc_id_r1, probs = seq(0, 1, 0.01), na.rm = TRUE)[99] # 98% percentile Query%

Target_90_perc = quantile(paralogues$hsapiens_paralog_perc_id, probs = seq(0, 1, 0.01), na.rm = TRUE)[91] # 90% percentile Target%

Query_90_perc = quantile(paralogues$hsapiens_paralog_perc_id_r1, probs = seq(0, 1, 0.01), na.rm = TRUE)[91] # 90% percentile Query%

Target_95_perc = quantile(paralogues$hsapiens_paralog_perc_id, probs = seq(0, 1, 0.01), na.rm = TRUE)[96] # 95% percentile Target%

Query_95_perc = quantile(paralogues$hsapiens_paralog_perc_id_r1, probs = seq(0, 1, 0.01), na.rm = TRUE)[96] # 95% percentile Query%

Target_99_perc = quantile(paralogues$hsapiens_paralog_perc_id, probs = seq(0, 1, 0.01), na.rm = TRUE)[100] # 99% percentile Target%

Query_99_perc = quantile(paralogues$hsapiens_paralog_perc_id_r1, probs = seq(0, 1, 0.01), na.rm = TRUE)[100] # 99% percentile Query%

Target_100_perc = quantile(paralogues$hsapiens_paralog_perc_id, probs = seq(0, 1, 0.01), na.rm = TRUE)[101] # 100% percentile Target%

Query_100_perc = quantile(paralogues$hsapiens_paralog_perc_id_r1, probs = seq(0, 1, 0.01), na.rm = TRUE)[101] # 100% percentile Query%


## Filter based on the percentiles

paralogues_98 = paralogues %>%
  filter(hsapiens_paralog_perc_id > Target_98_perc & hsapiens_paralog_perc_id_r1 > Query_98_perc) %>%
  rowwise() %>%
  mutate(paralogue = paste0(hsapiens_paralog_associated_gene_name, 
                            " (", as.character(round(hsapiens_paralog_perc_id)), ", ", 
                            as.character(round(hsapiens_paralog_perc_id_r1)), 
                            ", chr", hsapiens_paralog_chromosome, ")")) %>%
  group_by(ensembl_gene_id) %>%
  arrange(desc(hsapiens_paralog_perc_id)) %>%
  summarize(paralogues_98 = n())

paralogues_90 = paralogues %>%
  filter(hsapiens_paralog_perc_id > Target_90_perc & hsapiens_paralog_perc_id_r1 > Query_90_perc) %>%
  rowwise() %>%
  mutate(paralogue = paste0(hsapiens_paralog_associated_gene_name, 
                            " (", as.character(round(hsapiens_paralog_perc_id)), ", ", 
                            as.character(round(hsapiens_paralog_perc_id_r1)), 
                            ", chr", hsapiens_paralog_chromosome, ")")) %>%
  group_by(ensembl_gene_id) %>%
  arrange(desc(hsapiens_paralog_perc_id)) %>%
  summarize(paralogues_90 = n())

paralogues_95 = paralogues %>%
  filter(hsapiens_paralog_perc_id > Target_95_perc & hsapiens_paralog_perc_id_r1 > Query_95_perc) %>%
  rowwise() %>%
  mutate(paralogue = paste0(hsapiens_paralog_associated_gene_name, 
                            " (", as.character(round(hsapiens_paralog_perc_id)), ", ", 
                            as.character(round(hsapiens_paralog_perc_id_r1)), 
                            ", chr", hsapiens_paralog_chromosome, ")")) %>%
  group_by(ensembl_gene_id) %>%
  arrange(desc(hsapiens_paralog_perc_id)) %>%
  summarize(paralogues_95 = n())

paralogues_99 = paralogues %>%
  filter(hsapiens_paralog_perc_id > Target_99_perc & hsapiens_paralog_perc_id_r1 > Query_99_perc) %>%
  rowwise() %>%
  mutate(paralogue = paste0(hsapiens_paralog_associated_gene_name, 
                            " (", as.character(round(hsapiens_paralog_perc_id)), ", ", 
                            as.character(round(hsapiens_paralog_perc_id_r1)), 
                            ", chr", hsapiens_paralog_chromosome, ")")) %>%
  group_by(ensembl_gene_id) %>%
  arrange(desc(hsapiens_paralog_perc_id)) %>%
  summarize(paralogues_99 = n())

paralogues_100 = paralogues %>%
  filter(hsapiens_paralog_perc_id == Target_100_perc & hsapiens_paralog_perc_id_r1 == Query_100_perc) %>%
  rowwise() %>%
  mutate(paralogue = paste0(hsapiens_paralog_associated_gene_name, 
                            " (", as.character(round(hsapiens_paralog_perc_id)), ", ", 
                            as.character(round(hsapiens_paralog_perc_id_r1)), 
                            ", chr", hsapiens_paralog_chromosome, ")")) %>%
  group_by(ensembl_gene_id) %>%
  arrange(desc(hsapiens_paralog_perc_id)) %>%
  summarize(paralogues_100 = n())

# Data together

paralogues_final = paralogues_90 %>%
  left_join(paralogues_95) %>%
  left_join(paralogues_98) %>%
  left_join(paralogues_99) %>%
  left_join(paralogues_100) %>%
  replace(is.na(.), 0)


### Distance to telomere and centromere
############################################################################################################################

# downloaded the hg38 UCSC Table Browser the "gap" containing telomere coordinates
# downloaded the hg38 UCSC Table Browser the "centromeres" containing centromere coordinates

hg38_centromeres <- read_delim("raw/hg38_centromeres.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  group_by(chrom) %>%
  summarize(chrom = unique(chrom),
            start = min(chromStart),
            end = max(chromEnd))

hg38_gaps <- read_delim("raw/hg38_gaps.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  filter(type == "telomere")


## Gene positions 

gene_position = canonical_transcripts %>%
  dplyr::select(`Gene stable ID`, `HGNC ID`, `Chromosome/scaffold name`, `Transcription start site (TSS)`) %>%
  rename(chrom = "Chromosome/scaffold name",
         TSS = "Transcription start site (TSS)") %>%
  mutate(chrom = paste0("chr", chrom))


## Calculate distances

distances = list()
for (g in 1:nrow(gene_position)){
  chr = gene_position$chrom[g]
  pos = gene_position$TSS[g]
  
  centromere = hg38_centromeres %>% 
    filter(chrom == chr) %>%
    mutate(dist_centr = case_when(pos > end ~ pos - end,
                                  pos < start ~ start - pos),
           p_q = case_when(pos > end ~ "q",
                           pos < start ~ "p"))
  
  p_q = centromere$p_q
  telomeres = hg38_gaps %>%
    filter(chrom == chr) %>%
    mutate(dist_telom = case_when(p_q == "p" ~ pos - chromEnd,
                                  p_q == "q" ~ chromStart - pos)) %>%
    filter(dist_telom > 0) %>%
    pull(dist_telom)
  
  res = gene_position[g, ] %>%
    mutate(dist_centr = centromere$dist_centr, 
           dist_telom = telomeres)
  
  distances[[g]] = res
}

distances = do.call(rbind, distances)


### Dispensable genes
############################################################################################################################

# The separation of genes into classes was based on:
# - being a confirmed, PMT or no-disorder gene  
# - being related to a brain disorder (Neurologic_ClinicalSynopsis: TRUE or NA) 
# - the number of confident homozygous pLoF genotype observed (at least 1)  
# 
# The data for the last point was retrieved from Karczewski et al. 2020, Nature,581(7809):434-443
# Supplementary dataset 7 - List of genes where at least one confident homozygous pLoF genotype was observed
# Best current estimate of confidently LoF-tolerant genes based on the gnomAD dataset 


## Homozygous knockout genes

hom_ko_genes <- read_csv("raw/supplementary_dataset_7_hom_ko_genes.txt", 
                         col_names = FALSE)

hom_ko_genes_converted = convert_genes_to_approved(hom_ko_genes$X1)

hom_ko_genes_final = hom_ko_genes_converted[[1]]$approved_symbol

hom_ko_HGNC_OMIM = HGNC_OMIM %>% 
  mutate(hom_ko = case_when(`Approved symbol` %in% hom_ko_genes_final ~ 1, 
                            TRUE ~ 0)) %>%
  select(`Approved symbol`, monogenic_disorder, hom_ko, phen, phen_additional, Neurologic_ClinicalSynopsis)


## Gene classes

gene_classes = hom_ko_HGNC_OMIM %>%
  mutate(class = case_when(monogenic_disorder == "confirmed" & 
                             Neurologic_ClinicalSynopsis == TRUE & 
                             hom_ko == 0 ~ "training_validation_brain_disease", # = Cbi
                           monogenic_disorder == "no disorder" & 
                             is.na(Neurologic_ClinicalSynopsis) &
                             hom_ko == 1 ~ "training_validation_dispensable", # = NDt
                           
                           monogenic_disorder == "confirmed" & 
                             Neurologic_ClinicalSynopsis == TRUE & 
                             hom_ko == 1 ~ "confirmed_Brain_1homko", # = Cbt
                           monogenic_disorder == "confirmed" & 
                             is.na(Neurologic_ClinicalSynopsis) &
                             hom_ko == 0 ~ "confirmed_noBrain_0homko", # = C_i
                           monogenic_disorder == "confirmed" & 
                             is.na(Neurologic_ClinicalSynopsis) &
                             hom_ko == 1 ~ "confirmed_noBrain_1homko", # = C_t
                           
                           monogenic_disorder == "no disorder" & 
                             is.na(Neurologic_ClinicalSynopsis) &
                             hom_ko == 0 ~ "noDisorder_0homko", # = NDi
                           
                           monogenic_disorder == "PMT" & 
                             Neurologic_ClinicalSynopsis == TRUE & 
                             hom_ko == 0 ~ "PMT_Brain_0homko", # = PMTbi
                           monogenic_disorder == "PMT" & 
                             Neurologic_ClinicalSynopsis == TRUE & 
                             hom_ko == 1 ~ "PMT_Brain_1homko", # = PMTbt
                           monogenic_disorder == "PMT" &
                             is.na(Neurologic_ClinicalSynopsis) &
                             hom_ko == 0 ~ "PMT_noBrain_0homko", # = PMT_i
                           monogenic_disorder == "PMT" & 
                             is.na(Neurologic_ClinicalSynopsis) &
                             hom_ko == 1 ~ "PMT_noBrain_1homko")) # = PMT_t


### Final table
############################################################################################################################

ML_gnomad <- gnomad_preslice %>%
  select(ensembl_gene_id, obs_mis:max_af, num_coding_exons:gene_length) 

ML_GTEx <- expression_groups2 %>%
  rowwise() %>%
  mutate(ensemblID = str_split(Name, "[.]") %>% unlist() %>% .[1]) %>%
  select(-Name, -Description) %>%
  ungroup() %>%
  mutate(temp = rowMeans(.[, 1:(ncol(.)-1)], na.rm = TRUE)) %>%
  group_by(ensemblID) %>%
  arrange(desc(temp)) %>%
  slice(1) %>%
  select(-temp)

ML_tau <- tau_data %>%
  filter(!str_detect(Name, "PAR_Y")) %>%
  rowwise() %>%
  mutate(ensemblID = str_split(Name, "[.]") %>% unlist() %>% .[1]) %>%
  select(-Name, -Description) 

ML_Brainspan <- brainspan_data

ML_CpGdensity <- canonical_transcripts_2

ML_conserv_CDS_length <- conservation_scores_CDS_length %>% 
              select(-HGNC_symbol, -transcript_ID)

ML_paralogs <- paralogues_final
  
ML_distances <- distances %>%
  dplyr::select(`HGNC ID`, dist_centr, dist_telom)

ML_classes <- gene_classes %>%
  select(`Approved symbol`, class)


df = HGNC_OMIM %>%
  select(`HGNC ID`, `Approved symbol`, ensemblID, chr) %>%
  left_join(ML_classes) %>%
  select(-chr, chr) %>%
  filter(!is.na(chr)) %>% # remove the genes without chromosome or ensemblID assigned
  filter(!is.na(ensemblID)) %>%
  filter(`Approved symbol` != "MED14OS") %>% # in HGNC as protein coding but it as no coding exons !!!!!!!!!!!!!!!!
  left_join(ML_gnomad, by = c(ensemblID = "ensembl_gene_id")) %>%
  left_join(ML_conserv_CDS_length, c(ensemblID = "Gene stable ID")) %>%
  left_join(ML_CpGdensity, c(ensemblID = "Gene stable ID")) %>%
  left_join(ML_paralogs, by = c(ensemblID = "ensembl_gene_id")) %>%
  mutate(paralogues_90 = case_when(is.na(paralogues_90) ~ 0,
                                   TRUE ~ paralogues_90),
         paralogues_95 = case_when(is.na(paralogues_95) ~ 0,
                                   TRUE ~ paralogues_95),
         paralogues_98 = case_when(is.na(paralogues_98) ~ 0,
                                   TRUE ~ paralogues_98),
         paralogues_99 = case_when(is.na(paralogues_99) ~ 0,
                                   TRUE ~ paralogues_99),
         paralogues_100 = case_when(is.na(paralogues_100) ~ 0,
                                   TRUE ~ paralogues_100)) %>%
  left_join(ML_GTEx) %>%
  left_join(ML_tau) %>%
  left_join(ML_Brainspan, by = c(`Approved symbol` = "Gene")) %>%
  left_join(ML_distances) %>%
  unique() %>%
  select(-`HGNC ID`, -ensemblID) %>%
  mutate(class = case_when(class == "training_validation_brain_disease" ~ "Cbi",
                           class == "training_validation_dispensable" ~ "NDt",
                           class == "confirmed_Brain_1homko" ~ "Cbt",
                           class == "confirmed_noBrain_0homko" ~ "C_i",
                           class == "confirmed_noBrain_1homko" ~ "C_t",
                           class == "noDisorder_0homko" ~ "NDi",
                           class == "PMT_Brain_0homko" ~ "PMTbi",
                           class == "PMT_Brain_1homko" ~ "PMTbt",
                           class == "PMT_noBrain_0homko" ~ "PMT_i",
                           class == "PMT_noBrain_1homko" ~ "PMT_t"))

write_csv(df %>%
            select(-`HGNC ID`, -ensemblID), 
          "results/annotations_all_chr_ML/ML_data.csv", na = "")
