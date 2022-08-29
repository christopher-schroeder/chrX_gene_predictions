# author: "Elsa Leitao"

library(readr)
library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(seqinr)
library(biomaRt)
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
# GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
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
      print(gene)
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
  print(string)
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


### Annotation pseudoautosomal genes
############################################################################################################################

## Search done in HGNC (https://www.genenames.org/download/custom/)
# Curated by the HGNC: Approved symbol
# Select chromosomes: pseudoautosomal
# Select status: Approved
# https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&status=Approved&chr=XandY&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit
# saved as "raw/HGNC_PARgenes.txt"

HGNC_PARgenes <- read_delim("raw/HGNC_PARgenes.txt",
                            "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(pseudoautosomal = TRUE)


## Import HGNC genes with OMIM annotation from all chromosomes, restrict to chrX and annotate pseudoautosomal genes

HGNC_OMIM <- read_csv("results/HGNC_OMIM_all_chr/HGNCgenes_OMIMannotated.csv") %>%
  filter(Status == "Approved") %>%
  filter(chr == "X") %>%
  left_join(HGNC_PARgenes) %>%
  mutate(monogenic_disorder = case_when(gene_category %in% c("confirmed +CS", "confirmed -CS") ~ "confirmed",
                                        gene_category %in% c("PMT +CS", "PMT -CS") ~ "PMT",
                                        gene_category == "no disorder" ~ "no disorder")) %>%
  ungroup()


## Annotation Ensembl hg38 coordinates for genes in chrX 
############################################################################################################################

### BioMart (https://www.ensembl.org/biomart/martview/a210994a9cfc8d572ad222aa6e03331a)
## Database: Ensembl Genes 101
## Dataset: Human genes (GRCh38.p13)

## Filters:
# Chromosome/scaffold: X

## Attributes:
# Gene start (bp)
# Gene end (bp)
# Karyotype band
# Chromosome/scaffold name
# HGNC symbol
# HGNC ID
# Gene stable ID

# saved as "raw/mart_export_coordinates.txt"


ensembl <- read_delim("raw/mart_export_coordinates.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(`HGNC ID`, `Chromosome/scaffold name`, `Gene start (bp)`, `Gene end (bp)`, `Karyotype band`)


ensembl_coordinates = HGNC_OMIM %>%
  left_join(ensembl, by = c(`HGNC ID` = "HGNC ID")) %>%
  dplyr::select(`HGNC ID`, `Chromosome/scaffold name`, `Gene start (bp)`, `Gene end (bp)`, `Karyotype band`)


### gnomAD annotations for genes in chrX
############################################################################################################################

# file from "https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

# In OMIM tables, the gene names are HGNC approved symbols, but in gnomAD they aren't
# Before joining data, I converted gnomAD gene names to HGNC approved symbols. However, some of them had no approved symbols

# I kept only genes in chrX and the columns: approved_symbol, pLI, oe_lof_upper, syn_z, mis_z  

# In the gnomad table there where a few genes occuring more than once (multiple transcripts)
# I kept the occurrence with smaller oef_lof_upper (= LOEUF)


gnomad_lof_pre_slice = read.delim("raw/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz") %>%
  filter(chromosome == "X") %>%
  dplyr::select(gene, pLI, oe_lof_upper, syn_z, mis_z, exp_lof) %>%
  unique() %>%
  rowwise() %>%
  mutate(Gene = convert_single_gene_to_approved(gene)) %>%
  mutate(Gene = case_when(Gene == "" ~ as.character(gene),
                          TRUE ~ Gene)) %>%
  dplyr::select(-gene)

gnomad_lof = gnomad_lof_pre_slice %>%
  separate_rows(Gene, sep = ", ") %>%
  group_by(Gene) %>%
  arrange(oe_lof_upper) %>%
  dplyr::slice(1) %>%
  dplyr::rename(LOEUF = oe_lof_upper)


### GTEx annotations for genes in chrX
############################################################################################################################

# GTEx file "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz" 
# containing the median gene-level TPM by tissue (54 tissues) from all chromosomes

# (long calculation time (start))

GTEx_HGNC_approved <- read_delim("raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
                   "\t", escape_double = FALSE, trim_ws = TRUE,
                   skip = 2) %>% 
  rowwise() %>%
  mutate(Gene = convert_single_gene_to_approved(Description)) %>%
  mutate(Gene = case_when(Gene == "" ~ as.character(Description),
                          TRUE ~ Gene)) %>%
  select(-Description) %>%
  filter(!grepl("PAR_Y", Name)) # remove gene expression data from chrY in case of genes occurring in pseudoautosomal regions

# (long calculation time (end))


## Annotation of HGNC+OMIM with GTEX data

GTEX_annot_chrX_pre_sliced = HGNC_OMIM %>%
  left_join(GTEx_HGNC_approved, by = c(`Approved symbol` = "Gene")) %>%
  rowwise() %>%
  mutate(EnsemblID = str_split(Name, "[.]") %>% unlist() %>% .[1])


# 1 gene appears two times (2 different transcripts for the same gene)
# kept the row for which the two entrances of EnsemblID match


GTEX_annot_chrX = GTEX_annot_chrX_pre_sliced %>%
  mutate(ID_match = case_when(`Ensembl ID(supplied by Ensembl)` == EnsemblID ~ TRUE,
                              TRUE ~ FALSE)) %>%
  group_by(`Approved symbol`) %>%
  arrange(desc(ID_match)) %>%
  dplyr::slice(1)


### tau annotations for genes in chrX (using GTEx data)
############################################################################################################################

# "tau" as a measure of how tissue specific the expression of a gene is:
# - tau < 0.6 (broadly expressed)  
# - tau > 0.6 (expressed in a restricted number of tissues)  


## Calculate median of median expression in cerebellar tissues (2 regions) and in the other brain tissues (11 regions)
# Discard the expression in individual brain tissues

GTEx_chrX_tissues = GTEX_annot_chrX %>% 
  dplyr::select(`Approved symbol`, `Adipose - Subcutaneous`:`Whole Blood`) %>%
  rowwise() %>%
  mutate(Other_Brain_tissues = median(c(`Brain - Amygdala`,
                                        `Brain - Anterior cingulate cortex (BA24)`,
                                        `Brain - Caudate (basal ganglia)`, 
                                        `Brain - Cortex`,
                                        `Brain - Frontal Cortex (BA9)`,
                                        `Brain - Hippocampus`,
                                        `Brain - Hypothalamus`,
                                        `Brain - Nucleus accumbens (basal ganglia)`,
                                        `Brain - Putamen (basal ganglia)`,
                                        `Brain - Spinal cord (cervical c-1)`,
                                        `Brain - Substantia nigra`), 
                                      na.rm = TRUE),
         Cerebellar_tissues = median(c(`Brain - Cerebellar Hemisphere`,
                                       `Brain - Cerebellum`), 
                                     na.rm = TRUE)) %>%
  dplyr::select(-c(`Brain - Amygdala`,
                   `Brain - Anterior cingulate cortex (BA24)`,
                   `Brain - Caudate (basal ganglia)`,
                   `Brain - Cortex`,
                   `Brain - Frontal Cortex (BA9)`,
                   `Brain - Hippocampus`,
                   `Brain - Hypothalamus`,
                   `Brain - Nucleus accumbens (basal ganglia)`,
                   `Brain - Putamen (basal ganglia)`,
                   `Brain - Spinal cord (cervical c-1)`,
                   `Brain - Substantia nigra`,
                   `Brain - Cerebellar Hemisphere`,
                   `Brain - Cerebellum`)) 


## Calculate log2(TPM+1) median expression values for all tissues  
# Values shown in the columns "Other_Brain_tissues" and "Cerebellar_tissues" are log2(TPM+1)

GTEx_chrX_tissues_log <- GTEx_chrX_tissues %>%
  mutate_if(is.numeric, funs(log2(.+1)))


## Calculate the maximal expression value

GTEx_chrX_tissues_log = GTEx_chrX_tissues_log %>%
  gather(tissue, value, `Adipose - Subcutaneous`:`Cerebellar_tissues`) %>%
  group_by(`Approved symbol`) %>%
  summarize(max_expression = max(value)) %>%
  right_join(GTEx_chrX_tissues_log, by = "Approved symbol") %>%
  dplyr::select(-max_expression, max_expression) # put max_expression as last column


## Calculate expression profile component normalized by the maximal component value

GTEx_chrX_tissues_normlog = GTEx_chrX_tissues_log %>%
  mutate_if(is.numeric, funs(./max_expression)) %>%
  dplyr::select(-max_expression)


## Calculate tau

GTEx_chrX_tissues_normlog_1 = GTEx_chrX_tissues_normlog %>%
  mutate_if(is.numeric, funs(1-.)) 

GTEx_chrX_tissues_tau = GTEx_chrX_tissues_normlog_1 %>%
  gather(tissue, value, `Adipose - Subcutaneous`:`Cerebellar_tissues`) %>%
  group_by(`Approved symbol`) %>%
  summarize(tau = sum(value) / n()) %>%
  right_join(GTEx_chrX_tissues_normlog_1, by = "Approved symbol") %>%
  dplyr::select(`Approved symbol`, tau) %>%
  mutate(tissue_specificity = case_when(tau < 0.6 ~ "broadly expressed",
                                        TRUE ~ "tissue-specific"))


## Check the tissue having the maximum expression

GTEx_chrX_max_tissue = GTEx_chrX_tissues_normlog %>%
  left_join(GTEx_chrX_tissues_log %>% dplyr::select(`Approved symbol`, max_expression)) %>%
  filter(max_expression != 0) %>%
  gather(tissue, value, `Adipose - Subcutaneous`:`Cerebellar_tissues`) %>%
  group_by(`Approved symbol`) %>%
  arrange(desc(value)) %>%
  dplyr::slice(1) %>%
  summarize(max_tissue = tissue)


## Put GTEx data together
# When the highest expression was either in "Other_Brain_tissues" or "Cerebellar_tissues", 
# the genes was labeled TRUE in the column "brain_high"  

GTEx_chrX_allData = GTEx_chrX_tissues_tau %>%
  left_join(GTEx_chrX_tissues_log) %>%
  left_join(GTEx_chrX_max_tissue) %>%
  # left_join(GTEx_chrX_tissues_normlog)  %>%
  dplyr::select(`Approved symbol`, tau, max_expression, max_tissue, tissue_specificity, Other_Brain_tissues, Cerebellar_tissues) %>%
  dplyr::rename(max_expression_log2 = "max_expression") %>%
  mutate(Other_Brain_tissues = ifelse(is.nan(Other_Brain_tissues), NA, Other_Brain_tissues),
         Cerebellar_tissues = ifelse(is.nan(Cerebellar_tissues), NA, Cerebellar_tissues),
         tau = ifelse(is.nan(tau), NA, tau),
         brain_high = case_when(max_tissue %in% c("Other_Brain_tissues", "Cerebellar_tissues") ~ TRUE, 
                                TRUE ~ FALSE))


### X-inactivation annotations
############################################################################################################################

## Gather data from multiple publications related to genes escaping or not X-chromosome inactivation

# - Oliva et al. 2020 (https://science.sciencemag.org/content/369/6509/eaba3066.editor-summary) 
# - Tukiainen et al. 2017 (https://www.nature.com/articles/nature24265)  
# - Carrel and Willard 2005 (https://www.nature.com/articles/nature03479) (retrieved from Tukiainen and Zhang)  
# - Cotton et al. 2015 (https://academic.oup.com/hmg/article/24/6/1528/682766) (retrieved from Tukiainen)  
# - Wainer Katsir and Linial 2019 (https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5507-6)  
# - Zhang et al. 2013 (https://academic.oup.com/mbe/article/30/12/2588/1012600)  
# - Balaton et al. 2015 (https://bsd.biomedcentral.com/articles/10.1186/s13293-015-0053-7#Sec14) (summarizes 3 studies from Carrel and Cotton)  
# - Park et al. 2010 (https://academic.oup.com/mbe/article/27/11/2446/1120317) (based on Carrel genes; retrieved from Zhang)  
# - Schultz et al. 2015 (https://www.nature.com/articles/nature14465) (analysed the promoter of Carrel genes)  
#                        
# The tables retrieved from those publications are stored in folder "XCI":
# - Oliva_2020_aba3066-Table-S3.xlsx  
# - Tukiainen_2017_Suppl.Table.1.xlsx  
# - Tukiainen_2017_Suppl.Table.13.xlsx  
# - Wainer Katsir_Linia_2019_12864_2019_5507_MOESM6_ESM.xlsx  
# - Balaton_2015_13293_2015_53_MOESM1_ESM.xlsx  
# - Zhang_2013_2016_S3.xlsx  
#                        
# gene names converted to HGNC approved symbols  

# (long calculation time (start))

XCI_Oliva_2020_S3 <- read_excel("raw/XCI/Oliva_2020_aba3066-Table-S3.xlsx",
                                sheet = "GTEx v8 X-linked sex-biased gen",
                                na = "NA") %>%
  select(HUGO_gene_id, `Reported Escapee?`) %>%
  rename(Oliva2020_reported_escapee = "Reported Escapee?") %>%
  rowwise() %>%
  mutate(Gene = convert_single_gene_to_approved(HUGO_gene_id)) %>%
  mutate(Gene = case_when(Gene == "" ~ as.character(HUGO_gene_id),
                          TRUE ~ Gene)) %>%
  select(Gene, Oliva2020_reported_escapee) %>%
  unique()

Tukiainen_2017_S1 <- read_excel("raw/XCI/Tukiainen_2017_Suppl.Table.1.xlsx",
                                na = "NA", skip = 1) %>%
  select(`Gene name`, `Combined XCI status`, `XCI status...11`, `XCI status...14`) %>%
  rename(Tukiainen2017_previousCombinedXCIstatus = "Combined XCI status",
         Tukiainen2017_Cotton2015 = `XCI status...11`,
         Tukiainen2017_Carrel2005 = `XCI status...14`) %>%
  rowwise() %>%
  mutate(Gene = convert_single_gene_to_approved(`Gene name`)) %>%
  mutate(Gene = case_when(Gene == "" ~ as.character(`Gene name`),
                          TRUE ~ Gene)) %>%
  select(-`Gene name`) %>%
  unique()

Tukiainen_2017_S13 <- read_excel("raw/XCI/Tukiainen_2017_Suppl.Table.13.xlsx",
                                 na = "NA") %>%
  slice(1:(n()-1)) %>% # it has comments
  select(`Gene name`, `Reported XCI status`,
         `Sex-bias in GTEx`, `XCI across tissues`,
         `XCI in single cells`,
         `Possible new call for escape (support from several analyses)`) %>%
  rename(Tukiainen2017_ReportedXCIstatus = `Reported XCI status`,
         Tukiainen2017_SexBiasGTEx = `Sex-bias in GTEx`,
         Tukiainen2017_XCItissues = `XCI across tissues`,
         Tukiainen2017_XCIsingleCells = `XCI in single cells`,
         Tukiainen2017_NewCalls = `Possible new call for escape (support from several analyses)`) %>%
  rowwise() %>%
  mutate(Gene = convert_single_gene_to_approved(`Gene name`)) %>%
  mutate(Gene = case_when(Gene == "" ~ as.character(`Gene name`),
                          TRUE ~ Gene)) %>%
  select(-`Gene name`) %>%
  unique()

Wainer_Katsir_Linia_2019_S6 <- read_excel("raw/XCI/Wainer Katsir_Linia_2019_12864_2019_5507_MOESM6_ESM.xlsx",
                                          na = "NA", skip = 1) %>%
  select(`Gene Symbol`, Annotation) %>%
  rename(WainerKatsirLinia2019_Annotation = Annotation) %>%
  rowwise() %>%
  mutate(Gene = convert_single_gene_to_approved(`Gene Symbol`)) %>%
  mutate(Gene = case_when(Gene == "" ~ as.character(`Gene Symbol`),
                          TRUE ~ Gene)) %>%
  select(-`Gene Symbol`) %>%
  unique()

Balaton_2015_S1 <- read_excel("raw/XCI/Balaton_2015_13293_2015_53_MOESM1_ESM.xlsx",
                              sheet = "List") %>%
  select(`Gene Name`, `Balaton consensus calls`) %>%
  rename(Balaton2015_consensusCalls = "Balaton consensus calls") %>%
  rowwise() %>%
  mutate(Gene = convert_single_gene_to_approved(`Gene Name`)) %>%
  mutate(Gene = case_when(Gene == "" ~ as.character(`Gene Name`),
                          TRUE ~ Gene)) %>%
  select(-`Gene Name`) %>%
  unique()

Zhang_2013_2016_S3 <- read_excel("raw/XCI/Zhang_2013_2016_S3.xlsx",
                                 na = "NA") %>%
  select(Gene, `Our status`, `Merged Status`, `Carrel status`, `Park status`) %>%
  rename(gene = Gene,
         Zhang2013_Zhang = "Our status",
         Zhang2013_merged = "Merged Status",
         Zhang2013_Carrel = "Carrel status",
         Zhang2013_Park = "Park status") %>%
  rowwise() %>%
  mutate(Gene = convert_single_gene_to_approved(gene)) %>%
  mutate(Gene = case_when(Gene == "" ~ as.character(gene),
                          TRUE ~ Gene)) %>%
  select(-gene) %>%
  unique()

# (long calculation time (end))

## Gather all data into a single table 

XCI = XCI_Oliva_2020_S3 %>%
  full_join(Tukiainen_2017_S1, by = "Gene") %>%
  full_join(Tukiainen_2017_S13, by = "Gene") %>%
  full_join(Wainer_Katsir_Linia_2019_S6, by = "Gene") %>%
  full_join(Balaton_2015_S1, by = "Gene") %>%
  full_join(Zhang_2013_2016_S3, by = "Gene") %>%
  filter(Gene %in% approved_symbols$`Approved symbol`) %>%
  #to allow comparisons between columns
  mutate_at(vars(-("Gene")), str_to_lower) %>%
  mutate_all(str_replace_all, "inactivated", "inactive")
          

## Manual curation 

# 6 categories:  
# - high confidence escapee  
# - high confidence non-escapee  
# (all studies agreeing (perhaps except one study that could have done a manual mistake))  
#                        
# - low confidence escapee  
# - low confidence non-escapee  
# (disagreements but with a higher number of studies agreeing on one status)  
#                        
# - variable escapee  
# (most studies agreeing on variable escape)  
#                        
# - discordant  
# (similar number of studies agreeing on both status)  
#                        
# - na  
# (not enough data to have reliable evidence of XCI status)  
                       
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
# unzip to folder "raw/Brainspan_genes_matrix_csv"


## Import rows_metadata and convert gene names to approved symbols

# (long calculation time (start))

rows_metadata = read_csv("raw/Brainspan_genes_matrix_csv/rows_metadata.csv") %>%
  rowwise() %>%
  mutate(Gene = convert_single_gene_to_approved(gene_symbol)) %>% #converted the gene names to approved symbols
  mutate(Gene = case_when(Gene == "" ~ as.character(gene_symbol),
                          TRUE ~ Gene)) %>%
  select(-gene_symbol)

# (long calculation time (end))


## Correct a few conversions that were missed

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


## Keep data related to chrX genes (Row indexes referring to chrX genes) 

rows_chrX = which(rows_metadata$Gene %in% HGNC_OMIM$`Approved symbol`) 

rows_metadata_chrX = rows_metadata %>%
  filter(row_num %in% rows_chrX)

rows_metadata_chrX_curated = rows_metadata %>%
  filter(row_num %in% rows_chrX) %>%
  filter(ensembl_gene_id %in% HGNC_OMIM$`Ensembl gene ID` | entrez_id %in% HGNC_OMIM$`NCBI Gene ID(supplied by NCBI)`)%>%
  group_by(Gene) %>%
  arrange(desc(ensembl_gene_id %in% HGNC_OMIM$`Ensembl gene ID`)) %>%
  dplyr::slice(1)

genes_missing_5 = rows_metadata %>%
  filter(Gene %in% setdiff(rows_metadata_chrX$Gene %>% unique, rows_metadata_chrX_curated$Gene %>% unique()))

rows_metadata_chrX_curated2 = rows_metadata_chrX_curated %>%
  full_join(genes_missing_5)


## Expression data referring to chrX genes

expression_matrix_chrX = expression_matrix %>%
  filter(X1 %in% rows_metadata_chrX_curated2$row_num) %>%
  left_join(rows_metadata_chrX %>% 
              dplyr::select(Gene, row_num), by = c(X1 = "row_num"))


## Check columns_metadata to see which columns from the expression_matrix to summarize 

# Prenatal: age in pcw
prenatal_cols = paste0("X", (which(grepl("pcw", columns_metadata$age)) + 1)) 
# 1-237 (+ 1 because expression_matrix has a first column with the row_num)

# Post-natal 1: after birth until 4 year-old (inclusive)
postnatal_1_cols = paste0("X", which(columns_metadata$age %in% c("4 mos", "10 mos", "1 yrs",  "2 yrs", "3 yrs", "4 yrs")) + 1)
# 238-340 (+ 1 because expression_matrix has a first column with the row_num)

# Post-natal 2: older or equal to 8 years-old
postnatal_2_cols = paste0("X", which(columns_metadata$age %in% c("8 yrs", "11 yrs", "13 yrs", "15 yrs", "18 yrs", "19 yrs", "21 yrs", "23 yrs", "30 yrs",
                                                                 "36 yrs", "37 yrs", "40 yrs")) + 1)
# 341-524 (+ 1 because expression_matrix has a first column with the row_num)


## Simplify expression matrix by averaging on age groups while ignoring sex and tissue

expression_chrX_by_age = expression_matrix_chrX %>%
  dplyr::select(Gene, sample, TPM) %>%
  mutate(age_group = case_when(sample %in% prenatal_cols ~ "prenatal",
                               sample %in% postnatal_1_cols ~ "postnatal 1",
                               sample %in% postnatal_2_cols ~ "postnatal 2")) %>%
  mutate(age_group2 = case_when(sample %in% prenatal_cols ~ "prenatal",
                                sample %in% c(postnatal_1_cols, postnatal_2_cols) ~ "postnatal"))


## Median and mean

brainspan_median1 = expression_chrX_by_age %>%
  group_by(Gene, age_group) %>%
  summarize(median = median(TPM),
            mean = mean(TPM))

brainspan_median2 = expression_chrX_by_age %>%
  group_by(Gene, age_group2) %>%
  summarize(median = median(TPM),
            mean = mean(TPM)) %>%
  dplyr::rename(age_group = age_group2) %>%
  filter(age_group == "postnatal")

brainspan_median = brainspan_median1 %>% 
  rbind(brainspan_median2) %>%
  arrange(Gene)


## Mann-Whitney Prenatal vs Postnatal

MW1 = expression_chrX_by_age %>%
  dplyr::select(-sample) %>%
  group_by(Gene) %>%
  summarize(MW_Pre_Post = wilcox.test(TPM ~ age_group2)$p.value,
            MW_Pre_Post_adj = p.adjust(MW_Pre_Post, method = "bonferroni", n_groups(.)))


## Mann-Whitney Prenatal vs Postnatal 1

MW2 = expression_chrX_by_age %>%
  dplyr::select(-sample) %>%
  filter(age_group %in% c("prenatal", "postnatal 1")) %>%
  group_by(Gene) %>%
  summarize(MW_Pre_Post1 = wilcox.test(TPM ~ age_group)$p.value,
            MW_Pre_Post1_adj = p.adjust(MW_Pre_Post1, method = "bonferroni", n_groups(.)))


## Mann-Whitney Prenatal vs Postnatal 2

MW3 = expression_chrX_by_age %>%
  dplyr::select(-sample) %>%
  filter(age_group %in% c("prenatal", "postnatal 2")) %>%
  group_by(Gene) %>%
  summarize(MW_Pre_Post2 = wilcox.test(TPM ~ age_group)$p.value,
            MW_Pre_Post2_adj = p.adjust(MW_Pre_Post2, method = "bonferroni", n_groups(.)))


## Mann-Whitney Postnatal 1 vs Postnatal 2

MW4 = expression_chrX_by_age %>%
  dplyr::select(-sample) %>%
  filter(age_group %in% c("postnatal 1", "postnatal 2")) %>%
  group_by(Gene) %>%
  summarize(MW_Post1_Post2 = wilcox.test(TPM ~ age_group)$p.value,
            MW_Post1_Post2_adj = p.adjust(MW_Post1_Post2, method = "bonferroni", n_groups(.)))


## Combine data

brainspan_data <- brainspan_median %>%
  dplyr::select(-mean) %>%
  mutate(column = case_when(age_group == "prenatal" ~ "medianTPM_Pre",
                            age_group == "postnatal" ~ "medianTPM_Post",
                            age_group == "postnatal 1" ~ "medianTPM_Post1",
                            age_group == "postnatal 2" ~ "medianTPM_Post2")) %>%
  dplyr::select(-age_group) %>%
  spread(column, median) %>%
  dplyr::select(Gene, medianTPM_Pre, medianTPM_Post, medianTPM_Post1, medianTPM_Post2) %>%
  left_join(brainspan_median %>% 
              dplyr::select(-median) %>%
              mutate(column = case_when(age_group == "prenatal" ~ "meanTPM_Pre",
                                        age_group == "postnatal" ~ "meanTPM_Post",
                                        age_group == "postnatal 1" ~ "meanTPM_Post1",
                                        age_group == "postnatal 2" ~ "meanTPM_Post2")) %>%
              dplyr::select(-age_group) %>% 
              spread(column, mean) %>%
              dplyr::select(Gene, meanTPM_Pre, meanTPM_Post, meanTPM_Post1, meanTPM_Post2)) %>%
  left_join(MW1) %>%
  left_join(MW2) %>%
  left_join(MW3) %>%
  left_join(MW4) %>%
  mutate(median_logTPM_Pre = log((medianTPM_Pre + 1), 2),
         median_logTPM_Post = log((medianTPM_Post + 1), 2),
         median_logTPM_Post1 = log((medianTPM_Post1 + 1), 2),
         median_logTPM_Post2 = log((medianTPM_Post2 + 1), 2),
         mean_logTPM_Pre = log((meanTPM_Pre + 1), 2),
         mean_logTPM_Post = log((meanTPM_Post + 1), 2),
         mean_logTPM_Post1 = log((meanTPM_Post1 + 1), 2),
         mean_logTPM_Post2 = log((meanTPM_Post2 + 1), 2))


## Set categories

# Pre vs. Post:  
# - Pre_higher_vs_Post: adjusted p-value > 0.05 and higher mean TPM in Pre compared to Post  
# - Post_higher_vs_Pre: adjusted p-value > 0.05 and higher mean TPM in Post compared to Pre  
# - no_difference_Pre_Post  
# 
# Pre vs. Post1:  
# - Pre_higher_vs_Post1: adjusted p-value > 0.05 and higher mean TPM in Pre compared to Post1  
# - Post1_higher_vs_Pre: adjusted p-value > 0.05 and higher mean TPM in Post1 compared to Pre  
# - no_difference_Pre_Post1  
# 
# Pre vs. Post2:  
# - Pre_higher_vs_Post2: adjusted p-value > 0.05 and higher mean TPM in Pre compared to Post2  
# - Post2_higher_vs_Pre: adjusted p-value > 0.05 and higher mean TPM in Post2 compared to Pre  
# - no_difference_Pre_Post2  
# 
# Post1 vs. Post1:  
# - Post1_higher_vs_Post2: adjusted p-value > 0.05 and higher mean TPM in Post1 compared to Post2  
# - Post2_higher_vs_Post1: adjusted p-value > 0.05 and higher mean TPM in Post2 compared to Post1  
# - no_difference_Post1_Post2  

brainspan_data2 = brainspan_data %>% 
  mutate(BS_Pre_Post = case_when((MW_Pre_Post_adj < 0.05 & (meanTPM_Pre > 1 | meanTPM_Post > 1) & meanTPM_Pre > meanTPM_Post) ~ "Pre_higher_vs_Post",
                                 (MW_Pre_Post_adj < 0.05 & (meanTPM_Pre > 1 | meanTPM_Post > 1) & meanTPM_Post > meanTPM_Pre) ~ "Post_higher_vs_Pre",
                                 TRUE ~ "no_difference_Pre_Post"),
         BS_Pre_Post1 = case_when((MW_Pre_Post1_adj < 0.05 & (meanTPM_Pre > 1 | meanTPM_Post1 > 1) & meanTPM_Pre > meanTPM_Post1) ~ "Pre_higher_vs_Post1",
                                  (MW_Pre_Post1_adj < 0.05 & (meanTPM_Pre > 1 | meanTPM_Post1 > 1) & meanTPM_Post1 > meanTPM_Pre) ~ "Post1_higher_vs_Pre",
                                  TRUE ~ "no difference_Pre_Post1"),
         BS_Pre_Post2 = case_when((MW_Pre_Post2_adj < 0.05 & (meanTPM_Pre > 1 | meanTPM_Post2 > 1) & meanTPM_Pre > meanTPM_Post2) ~ "Pre_higher_vs_Post2",
                                  (MW_Pre_Post2_adj < 0.05 & (meanTPM_Pre > 1 | meanTPM_Post2 > 1) & meanTPM_Post2 > meanTPM_Pre) ~ "Post2_higher_vs_Pre",
                                  TRUE ~ "no difference_Pre_Post2"),
         BS_Post1_Post2 = case_when((MW_Post1_Post2_adj < 0.05 & (meanTPM_Post1 > 1 | meanTPM_Post2 > 1) & meanTPM_Post1 > meanTPM_Post2) ~ "Post1_higher_vs_Post2",
                                    (MW_Post1_Post2_adj < 0.05 & (meanTPM_Post1 > 1 | meanTPM_Post2 > 1) & meanTPM_Post2 > meanTPM_Post1) ~ "Post2_higher_vs_Post1",
                                    TRUE ~ "no difference_Post1_Post2")) 

  
### Uniprot annotations
############################################################################################################################

## Search done in Uniprot:
# https://www.uniprot.org/uniprot/?query=&fil=reviewed%3Ayes%20AND%20organism%3A%22Homo%20sapiens%20(Human)%20%5B9606%5D%22&columns=id%2Centry%20name%2Creviewed%2Ccomment(FUNCTION)%2Cgo(biological%20process)%2Cgo(molecular%20function)%2Ccomment(SUBCELLULAR%20LOCATION)%2Ccomment(SUBUNIT%20STRUCTURE)
# - Entry name  
# - Function [CC]  
# - Subunit structure [CC]  
# - Gene ontology (biological process)  
# - Gene ontology (molecular function)  
# - Subcellular location [CC]  
# saved as "uniprot_reviewed_human.xlsx"


uniprot <- read_excel("raw/uniprot_reviewed_human.xlsx") %>%
  dplyr::select(-`Entry name`, -Status)


uniprot_chrX = HGNC_OMIM %>%
  dplyr::select(`Approved symbol`, `UniProt ID(supplied by UniProt)`) %>%
  separate_rows(`UniProt ID(supplied by UniProt)`, sep = ", ") %>%
  left_join(uniprot, by = c(`UniProt ID(supplied by UniProt)` = "Entry")) %>%
  group_by(`Approved symbol`) %>%
  summarize(`UniProt ID(supplied by UniProt)` = paste0(`UniProt ID(supplied by UniProt)`[!is.na(`UniProt ID(supplied by UniProt)`)], collapse = ", "),
            `Function [CC]` = paste0(`Function [CC]`[!is.na(`Function [CC]`)], collapse = "//"),
            `Gene ontology (biological process)` = paste0(`Gene ontology (biological process)`[!is.na(`Gene ontology (biological process)`)], collapse = "//"),
            `Gene ontology (molecular function)` = paste0(`Gene ontology (molecular function)`[!is.na(`Gene ontology (molecular function)`)], collapse = "//"),
            `Subcellular location [CC]` = paste0(`Subcellular location [CC]`[!is.na(`Subcellular location [CC]`)], collapse = "//"),
            `Subunit structure [CC]` = paste0(`Subunit structure [CC]`[!is.na(`Subunit structure [CC]`)], collapse = "//"))


### Promoter CpG density
############################################################################################################################

## Transcripts retrieved from chrX from Ensembl via BioMart using the filters/attributes below:
  
# Ensembl Genes 102  
# Dataset:  
# - Human genes (GRCh38.p13)  
# Filters:  
# - Chromosome/scaffold: X  
# Attributes:  
# - Gene stable ID  
# - Transcript stable ID  
# - Transcript stable ID version  
# - HGNC ID  
# - HGNC symbol  
# - Gene start (bp)  
# - Gene end (bp)  
# - Transcript start (bp)  
# - Transcript end (bp)  
# - Transcription start site (TSS)  
# - GENCODE basic annotation  
# - APPRIS annotation  
# - RefSeq match transcript  
# saved as "raw/mart_export_transcripts_all_2.txt"  

# For restricting to canonical transcripts I considered MANE and APPRIS annotations:  
#   - MANE Select transcripts (Matched Annotation between NCBI and EBI) were independently identified by both Ensembl and NCBI as the most biologically relevant.   
# - APPRIS is a system to annotate alternatively spliced transcripts based on a range of computational methods.  
# 
# More information about the transcript quality tags can be found in https://www.ensembl.org/info/genome/genebuild/transcript_quality_tags.html.

# Note: for one case I had to use the "Ensembl ID(supplied by Ensembl)" and not the HGNC curated version because it was wrong.  
# 
# The ratio observed/expected (Obs/Exp or O/E) CpG was calculated as follows: Obs/Exp CpG = (Number of CpG / (Number of C x Number of G)) x N, where N is the total number of nucleotides in the sequence being analysed.  

ensembl_transcripts <- read_delim("raw/mart_export_transcripts_all_2.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(length = `Transcript end (bp)` - `Transcript start (bp)`)


## MANE and APPRIS annotations
# - MANE Select transcripts (Matched Annotation between NCBI and EBI) were independently identified by both Ensembl and NCBI as the most biologically relevant.   
# - APPRIS is a system to annotate alternatively spliced transcripts based on a range of computational methods.
MANE_transcripts = ensembl_transcripts %>% 
  filter(!is.na(`RefSeq match transcript`) & `Gene stable ID` %in% HGNC_OMIM$`Ensembl gene ID`)

APPRIS_transcripts = ensembl_transcripts %>% 
  filter(!`Gene stable ID` %in% MANE_transcripts$`Gene stable ID` & !is.na(`APPRIS annotation`) & `Gene stable ID` %in% HGNC_OMIM$`Ensembl gene ID`) %>%
  group_by(`Gene stable ID`) %>%
  arrange(desc(length)) %>%
  dplyr::slice(1)

other_transcripts = ensembl_transcripts %>%
  filter(`Gene stable ID` %in% HGNC_OMIM$`Ensembl gene ID` &
           !`Gene stable ID` %in% MANE_transcripts$`Gene stable ID` &
           !`Gene stable ID` %in% APPRIS_transcripts$`Gene stable ID`)

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
  dplyr::select(`Transcription start site (TSS)`) %>%
  mutate(chr = "chrX",
         start = `Transcription start site (TSS)`- 2000,
         end = `Transcription start site (TSS)`+ 2000) %>%
  dplyr::select(chr, start, end) %>%
  write_tsv( "results/annotations_chrX/promoters_chrX.bed", col_names = FALSE)


## Get promoter sequences
# "results/annotations_chrX/promoters_chrX.bed" uploaded as a "Custom Track" in UCSC and obtained the sequences through "Table Browser"
# saved as "raw/promoter_chrX_seqs.txt"

promoter_seqs = read.fasta("raw/promoter_chrX_seqs.txt", as.string = TRUE)


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
  mutate(chr = "chrX",
         start = `Transcription start site (TSS)`- 1999,
         end = `Transcription start site (TSS)`+ 2000) %>%
  unite(tmp1, chr:start, sep = ":") %>%
  unite(promoter, tmp1:end, sep = "-") %>%
  left_join(CpG_density_list) %>%
  dplyr::select(`HGNC symbol`, CpG_density)


### Promoter and exon-conservation across species
############################################################################################################################

## Use conservation scores from phastCons100way.UCSC.hg38
# PhastCons score ranges from 0 to 1 and represents the probability that a given nucleotide is conserved

phast <- phastCons100way.UCSC.hg38


## Exons

# biomaRt 
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl", 
                   host = "https://www.ensembl.org")


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
  filter(!is.na(cds_start) & !is.na(cds_end)) %>% # removes non-coding exons
  mutate(start = case_when(all_coding == FALSE & cds_start == 1 & exon_chrom_start == min(.$exon_chrom_start) ~ exon_chrom_end - cds_end + cds_start,  #correct for partially coding/non-coding exons
                           all_coding == FALSE & cds_start != 1 & exon_chrom_start == min(.$exon_chrom_start) ~ exon_chrom_end - cds_end + cds_start,
                           TRUE ~ exon_chrom_start),
         end = case_when(all_coding == FALSE & cds_start != 1 & exon_chrom_end == max(.$exon_chrom_end) ~ exon_chrom_start + cds_end - cds_start,
                         all_coding == FALSE & cds_start == 1 & exon_chrom_end == max(.$exon_chrom_end) ~ exon_chrom_start + cds_end - cds_start,
                         TRUE ~ exon_chrom_end))


# Scores from coding exons

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
  dplyr::select(`Transcript stable ID`, `HGNC symbol`, `Transcription start site (TSS)`) %>%
  mutate(chr = "chrX",
         start = `Transcription start site (TSS)`- 2000,
         end = `Transcription start site (TSS)`+ 2000)

promoter_scores = list()
for(i in 1:nrow(promoter_coordinates)){
  chr = promoter_coordinates$chr[i]
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

# ((long calculation time (end)))


## Both conservation scores

conservation_scores = exon_scores_2 %>%
  left_join(promoter_scores_2)


### CDS_length
############################################################################################################################

length = exons_coordinates %>%
  mutate(CDS_length = end - start + 1) %>%
  group_by(`HGNC symbol`) %>%
  summarize(CDS_length = sum(CDS_length))


## Data together

conservation_scores_CDS_length = conservation_scores %>%
  left_join(length, by = c(HGNC_symbol = "HGNC symbol"))


### Paralogues from Ensembl
############################################################################################################################

# "Paralogues are defined in Ensembl as genes for which the most common ancestor node is a duplication event
# These ancestral duplications are represented by red nodes in the gene trees.
# The table shows the taxonomic level of the ancestor duplication node,
# the Ensembl gene ID and name, the location of the paralogue, 
# and the percent of identical amino acids in the paralogue compared with the gene of interest (Target %ID)
# The identity of the gene of interest when compared with the paralogue is the query %ID."


## Bring paralogues genes on chrX using biomaRt  

paralogues = getBM(attributes = c("ensembl_gene_id",
                                  "hsapiens_paralog_associated_gene_name",
                                  "hsapiens_paralog_chromosome",
                                  "hsapiens_paralog_orthology_type",
                                  "hsapiens_paralog_subtype",
                                  "hsapiens_paralog_perc_id",
                                  "hsapiens_paralog_perc_id_r1"),
                   filters = c("ensembl_gene_id"),
                   values = list(canonical_transcripts$`Gene stable ID`),
                   mart = ensembl) %>%
  left_join(canonical_transcripts %>%
              dplyr::select(`HGNC symbol`, `Gene stable ID`), by = c(ensembl_gene_id = "Gene stable ID"))


## Calculate the 90 and 95 percentiles of the Target% and Query%

Target_90_perc = quantile(paralogues$hsapiens_paralog_perc_id, probs = seq(0, 1, 0.05), na.rm = TRUE)[19] # 90% percentile Target%

Target_95_perc = quantile(paralogues$hsapiens_paralog_perc_id, probs = seq(0, 1, 0.05), na.rm = TRUE)[20] # 95% percentile Target%

Query_90_perc = quantile(paralogues$hsapiens_paralog_perc_id_r1, probs = seq(0, 1, 0.05), na.rm = TRUE)[19] # 90% percentile Query%

Query_95_perc = quantile(paralogues$hsapiens_paralog_perc_id_r1, probs = seq(0, 1, 0.05), na.rm = TRUE)[20] # 95% percentile Query%


## Filter based on the percentiles

paralogues_90 = paralogues %>%
  filter(hsapiens_paralog_perc_id > Target_90_perc & hsapiens_paralog_perc_id_r1 > Query_90_perc) %>%
  rowwise() %>%
  mutate(paralogue = paste0(hsapiens_paralog_associated_gene_name, 
                            " (", as.character(round(hsapiens_paralog_perc_id)), ", ", 
                            as.character(round(hsapiens_paralog_perc_id_r1)), 
                            ", chr", hsapiens_paralog_chromosome, ")")) %>%
  group_by(ensembl_gene_id) %>%
  arrange(desc(hsapiens_paralog_perc_id)) %>%
  summarize(number_paralogues_90_perc = n(),
            paralogues_90_perc = paste0(paralogue, collapse = "; "),
            paralogues_90_logic = TRUE) 

paralogues_95 = paralogues %>%
  filter(hsapiens_paralog_perc_id > Target_95_perc & hsapiens_paralog_perc_id_r1 > Query_95_perc) %>%
  rowwise() %>%
  mutate(paralogue = paste0(hsapiens_paralog_associated_gene_name, 
                            " (", as.character(round(hsapiens_paralog_perc_id)), ", ", 
                            as.character(round(hsapiens_paralog_perc_id_r1)), 
                            ", chr", hsapiens_paralog_chromosome, ")")) %>%
  group_by(ensembl_gene_id) %>%
  arrange(desc(hsapiens_paralog_perc_id)) %>%
  summarize(number_paralogues_95_perc = n(),
            paralogues_95_perc = paste0(paralogue, collapse = "; "),
            paralogues_95_logic = TRUE)


### Final table
############################################################################################################################

df = HGNC_OMIM %>%
  filter(`Approved symbol` != "MED14OS") %>% # in HGNC as protein coding but it as no coding exons 
  left_join(uniprot_chrX, by = c("Approved symbol", "UniProt ID(supplied by UniProt)")) %>%
  left_join(ensembl_coordinates, by = c("HGNC ID")) %>%
  left_join(gnomad_lof, by = c(`Approved symbol` = "Gene")) %>%
  left_join(canonical_transcripts_2, by = c(`Approved symbol` = "HGNC symbol")) %>%
  left_join(conservation_scores_CDS_length, c(`Approved symbol` = "HGNC_symbol")) %>%
  left_join(paralogues_90, c(`Ensembl gene ID` = "ensembl_gene_id")) %>%
  left_join(paralogues_95, c(`Ensembl gene ID` = "ensembl_gene_id")) %>%
  left_join(GTEx_chrX_allData) %>%
  left_join(XCI_curated_2, by = c(`Approved symbol` = "Gene")) %>%
  left_join(brainspan_data2 %>%
              dplyr::select(Gene, BS_Pre_Post, BS_Pre_Post1, BS_Pre_Post2, BS_Post1_Post2), by = c(`Approved symbol` = "Gene")) %>%
  mutate(ClinicalSynopsis = case_when(ClinicalSynopsis == T ~ TRUE,
                                      TRUE ~ FALSE))%>%
  mutate(Neurologic_ClinicalSynopsis = case_when(Neurologic_ClinicalSynopsis == T ~ TRUE,
                                                 TRUE ~ FALSE)) %>%
  mutate(intellectual_disability = case_when(intellectual_disability == T ~ TRUE,
                                             TRUE ~ FALSE)) %>%
  mutate(seizures = case_when(seizures == T ~ TRUE,
                              TRUE ~ FALSE)) %>%
  mutate(motor_development = case_when(motor_development == T ~ TRUE,
                                       TRUE ~ FALSE)) %>%
  mutate(language_development = case_when(language_development == T ~ TRUE,
                                          TRUE ~ FALSE)) %>%
  mutate(spasticity = case_when(spasticity == T ~ TRUE,
                                TRUE ~ FALSE)) %>%
  mutate(ataxia = case_when(ataxia == T ~ TRUE,
                            TRUE ~ FALSE)) %>%
  mutate(pseudoautosomal = case_when(pseudoautosomal == T ~ TRUE,
                                     TRUE ~ FALSE)) %>%
  mutate(gene_group = case_when(monogenic_disorder == "confirmed" ~ "confirmed monogenic disorder",
                                monogenic_disorder == "PMT" ~ "PMT",
                                TRUE ~ "no disorder"),
         gene_group = factor(gene_group, levels = c("confirmed monogenic disorder",
                                                    "PMT",
                                                    "no disorder")))

write_csv(df, "results/annotations_chrX/HGNC_OMIM_annotated_plus_multiple_annotations.csv", na = "")


# NOTE: MED14OS is in HGNC as protein-coding but it has no coding exons in Ensembl (it is a lncRNA)
# It has been curated in Uniprot as "Product of a dubious CDS prediction"
# This gene was excluded from the analysis