# author: "Elsa Leitao"

library(UniprotR)
library(dplyr)
library(stringr)
library(biomaRt)
library(ggplot2)
library(tidyr)
library(ggforce)
library(readxl)
library(grid)


### Functions
############################################################################################################################

## Convert vector with domains and regions from uniprot into a df

convert_uniprot_vector_df = function(vector_uniprot, type){
  df = list()
  for(i in 1:length(vector_uniprot)){
    split = vector_uniprot[i] %>%
      str_split("[/]") %>%
      unlist()
    
    start = split[1] %>%
      str_split("[..]") %>%
      unlist() %>%
      .[1] %>%
      as.integer()
    
    end = split[1] %>%
      str_split("[..]") %>%
      unlist() %>%
      .[3] %>%
      trimws() %>%
      str_replace("[;]", "") %>%
      as.integer()
    
    name = split[2] %>%
      str_replace("note=", "") %>%
      trimws() %>%
      str_replace("[;]", "")
    
    df[[i]] = data.frame(protein = i, 
                         name = name,
                         start = start,
                         end = end,
                         type = type)
  }
  
  df = do.call(rbind, df)
  
  return(df)
}


### Import information for gene schemes
############################################################################################################################

# Selected genes

gene_for_scheme <- read_excel("raw/gene_for_scheme.xlsx") %>%
  unique() 


# Selected genes validation

gene_for_scheme_val <- read_excel("raw/gene_for_scheme_validation.xlsx") %>%
  unique() 


### Bring exon sizes from Ensembl
############################################################################################################################

## RefSeq transcript ID

# Selected genes

refIDS = gene_for_scheme %>%
  dplyr::select(RefSeq_transcript_id) %>%
  rowwise() %>%
  mutate(id = str_split(RefSeq_transcript_id, "[.]") %>% unlist() %>% .[1])


# Selected genes validation

refIDS_val = gene_for_scheme_val %>%
  dplyr::select(RefSeq_transcript_id) %>%
  rowwise() %>%
  mutate(id = str_split(RefSeq_transcript_id, "[.]") %>% unlist() %>% .[1])


## Bring coding exons coordinates from selected transcripts

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")

# Selected genes

exons_coordinates = getBM(attributes = c("external_gene_name", "ensembl_transcript_id",
                                         "cds_start",
                                         "cds_end"),
                          filters = c("refseq_mrna"),
                          values = list(refIDS %>% pull(id) %>% unique()),
                          mart = ensembl) %>%
  group_by(ensembl_transcript_id) %>%
  arrange(cds_start) %>%
  mutate(name = 1:n(),
         type = "exon") %>%
  rename(start = "cds_start",
         end = "cds_end") %>%
  na.omit() # due to exons that are non-coding


# Selected genes validation

exons_coordinates_val = getBM(attributes = c("external_gene_name", "ensembl_transcript_id",
                                             "cds_start",
                                             "cds_end"),
                              filters = c("refseq_mrna"),
                              values = list(refIDS_val %>% pull(id) %>% unique()),
                              mart = ensembl) %>%
  group_by(ensembl_transcript_id) %>%
  arrange(cds_start) %>%
  mutate(name = 1:n(),
         type = "exon") %>%
  rename(start = "cds_start",
         end = "cds_end") %>%
  na.omit() # due to exons that are non-coding


### Uniprot data
############################################################################################################################

# Selected genes

domains_list = list()
zn_fingers_list = list()
repeats_list = list()
coiledCoil_list = list()
tm_list = list()

for (p in gene_for_scheme$Uniprot_id %>% unique()) {
  
  # Protein domains, Zinc fingers, Repeats, coiled coils 
  uniprot = GetFamily_Domains(p, directorypath = NULL)
  
  domains = uniprot %>%
    pull(Domain..FT.) %>%
    str_split("DOMAIN ") %>%
    unlist() %>%
    .[-1]
  
  zn_fingers = uniprot %>%
    pull(Zinc.finger) %>%
    str_split("ZN_FING ") %>%
    unlist() %>%
    .[-1]
  
  repeats = uniprot %>%
    pull(Repeat) %>%
    str_split("REPEAT ") %>%
    unlist() %>%
    .[-1] %>%
    str_subset("NHL|YD|WD")
  
  coiledCoil = uniprot %>%
    pull(Coiled.coil) %>%
    str_split("COILED ") %>%
    unlist() %>%
    .[-1] 
  
  
  #Transmembrane domain
  uniprot2 = GetSubcellular_location(p, directorypath = NULL)
  
  tm = uniprot2 %>%
    pull(Transmembrane) %>%
    str_split("TRANSMEM ") %>%
    unlist() %>%
    .[-1]
  
  
  # Get info
  
  if(length(domains) == 0){
    df_domains = data.frame(protein = p,
                            name = NA_character_,
                            start = NA_integer_,
                            end = NA_integer_,
                            type = "none")
  } else {
    df_domains = convert_uniprot_vector_df(domains, "domain") %>%
      mutate(protein = p)
  }
  
  if(length(zn_fingers) == 0){
    df_zn_fingers = data.frame(protein = p,
                               name = NA_character_,
                               start = NA_integer_,
                               end = NA_integer_,
                               type = "none")
  } else {
    df_zn_fingers = convert_uniprot_vector_df(zn_fingers, "zn_fingers") %>%
      mutate(protein = p)
  }
  
  if(length(repeats) == 0){
    df_repeats = data.frame(protein = p,
                            name = NA_character_,
                            start = NA_integer_,
                            end = NA_integer_,
                            type = "none")
  } else {
    df_repeats = convert_uniprot_vector_df(repeats, "NHL_YD_WD_repeats") %>%
      mutate(protein = p)
  }
  
  if(length(coiledCoil) == 0){
    df_coiledCoil = data.frame(protein = p,
                               name = NA_character_,
                               start = NA_integer_,
                               end = NA_integer_,
                               type = "none")
  } else {
    df_coiledCoil = convert_uniprot_vector_df(coiledCoil, "coiledCoil") %>%
      mutate(protein = p)
  }
  
  if(length(tm) == 0){
    df_tm = data.frame(protein = p,
                       name = NA_character_,
                       start = NA_integer_,
                       end = NA_integer_,
                       type = "none")
  } else {
    df_tm = convert_uniprot_vector_df(tm, "transmembrane") %>%
      mutate(protein = p)
  }
  domains_list [[p]]= df_domains
  zn_fingers_list [[p]]= df_zn_fingers
  repeats_list [[p]]= df_repeats
  coiledCoil_list [[p]]= df_coiledCoil
  tm_list [[p]]= df_tm
}

domains_list = do.call(rbind, domains_list)
zn_fingers_list = do.call(rbind, zn_fingers_list)
repeats_list = do.call(rbind, repeats_list)
coiledCoil_list = do.call(rbind, coiledCoil_list)
tm_list = do.call(rbind, tm_list)


# Selected genes validation

domains_list_val = list()
repeats_list_val = list()
tm_list_val = list()

for (p in gene_for_scheme_val$Uniprot_id %>% unique()) {
  
  uniprot = GetFamily_Domains(p, directorypath = NULL)
  
  domains = uniprot %>%
    pull(Domain..FT.) %>%
    str_split("DOMAIN ") %>%
    unlist() %>%
    .[-1]
  
  repeats = uniprot %>%
    pull(Repeat) %>%
    str_split("REPEAT ") %>%
    unlist() %>%
    .[-1] %>%
    str_subset("ANK")
  
  #Transmembrane domain
  uniprot2 = GetSubcellular_location(p, directorypath = NULL)
  
  tm = uniprot2 %>%
    pull(Transmembrane) %>%
    str_split("TRANSMEM ") %>%
    unlist() %>%
    .[-1]
  

  if(length(domains) == 0){
    df_domains = data.frame(protein = p,
                            name = NA_character_,
                            start = NA_integer_,
                            end = NA_integer_,
                            type = "none")
  } else {
    df_domains = convert_uniprot_vector_df(domains, "domain") %>%
      mutate(protein = p)
  }
  
  if(length(repeats) == 0){
    df_repeats = data.frame(protein = p,
                            name = NA_character_,
                            start = NA_integer_,
                            end = NA_integer_,
                            type = "none")
  } else {
    df_repeats = convert_uniprot_vector_df(repeats, "ANK_repeats") %>%
      mutate(protein = p)
  }
  
  if(length(tm) == 0){
    df_tm = data.frame(protein = p,
                       name = NA_character_,
                       start = NA_integer_,
                       end = NA_integer_,
                       type = "none")
  } else {
    df_tm = convert_uniprot_vector_df(tm, "transmembrane") %>%
      mutate(protein = p)
  }
  domains_list_val [[p]]= df_domains
  repeats_list_val [[p]]= df_repeats
  tm_list_val [[p]]= df_tm
  
}

domains_list_val = do.call(rbind, domains_list_val)
repeats_list_val = do.call(rbind, repeats_list_val)
tm_list_val = do.call(rbind, tm_list_val)


### Gene schemes with variants
############################################################################################################################

## Color codes

colors = c(adjustcolor("#ff0000", alpha.f = 0.2), 
           adjustcolor("#ff9900", alpha.f = 0.2), 
           adjustcolor("#00b04f", alpha.f = 0.2), 
           adjustcolor("#0070c0", alpha.f = 0.2), 
           adjustcolor("#6f30a0", alpha.f = 0.2), 
           adjustcolor("#c00000", alpha.f = 0.2),
           adjustcolor("#ffbf00", alpha.f = 0.2),
           adjustcolor("#92d050", alpha.f = 0.2),
           adjustcolor("#00b0f0", alpha.f = 0.2),
           adjustcolor("#002060", alpha.f = 0.2))

colors_regions = c(adjustcolor("#ff0000", alpha.f = 0.8),
                   adjustcolor("#ffff00", alpha.f = 0.8),
                   adjustcolor("#00b04f", alpha.f = 0.8),
                   adjustcolor("#0070c0", alpha.f = 0.8),
                   adjustcolor("#6f30a0", alpha.f = 0.8),
                   adjustcolor("#c00000", alpha.f = 0.8),
                   adjustcolor("#ffbf00", alpha.f = 0.8),
                   adjustcolor("#92d050", alpha.f = 0.8),
                   adjustcolor("#00b0f0", alpha.f = 0.8),
                   adjustcolor("#002060", alpha.f = 0.8))


## Selected genes TENM1 (half the scale compared to TENM1)

for (g in unique(gene_for_scheme$gene_name) %>% .[! . %in% c("TENM1")]){
  print(g)
  
  tbl = gene_for_scheme %>%
    filter(gene_name == g)
  
  transcript = tbl$RefSeq_transcript_id %>% 
    unique()
  
  protein_id = tbl$Uniprot_id %>%
    unique()
  
  exons = exons_coordinates %>%
    filter(external_gene_name == g) %>%
    mutate(name = as.character(name)) %>%
    dplyr::select(name, start, end, type)
  
  domains = domains_list %>%
    filter(protein == protein_id) %>%
    dplyr::select(name, start, end, type)
  
  zn_fingers = zn_fingers_list %>%
    filter(protein == protein_id) %>%
    dplyr::select(name, start, end, type)
  
  repeats = repeats_list %>%
    filter(protein == protein_id) %>%
    dplyr::select(name, start, end, type)
  
  coiledCoil = coiledCoil_list %>%
    filter(protein == protein_id) %>%
    dplyr::select(name, start, end, type)
  
  tm = tm_list %>%
    filter(protein == protein_id) %>%
    dplyr::select(name, start, end, type)
  
  cds_seq_size = max(exons$end) - min(exons$start) + 1 - 3 # without the stop codon
  protein_size = cds_seq_size / 3
  
  fragments = rbind(exons, domains, zn_fingers, repeats, coiledCoil, tm) %>%
    mutate(position = end)
  
  position_max = (max(exons_coordinates$end) + 200) / 2
  
  graph_gene = fragments %>%
    ggplot(aes(position)) +
    
    ### Add transcript and protein IDs
    labs(title = substitute(italic(x)/y/z, list(x = g, y = transcript, z = protein_id))) +
    coord_cartesian(clip = "off") +
    xlim(0, position_max) + ###
    theme_classic() +
    theme(
      plot.title = element_text(size = 20),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    
    ### Add transcript rectangle and first and last positions
    annotate("rect", xmin = 1, xmax = cds_seq_size, ymin = cds_seq_size/100*1.5, ymax = cds_seq_size/100*2.5,
             color = "black", fill = "white") +
    annotate("text", x = cds_seq_size + cds_seq_size*0.02, y = cds_seq_size/100*2, label = paste0(cds_seq_size, " bp"),
             hjust = 0) +
    
    ### Add protein rectangle and first and last positions
    annotate("rect", xmin = 1, xmax = cds_seq_size, ymin = 0, ymax = cds_seq_size/100,
             color = "black", fill = "white") +
    annotate("text", x = cds_seq_size + cds_seq_size*0.02, y = cds_seq_size/100/2, label = paste0(protein_size, " aa"),
             hjust = 0)
  
  
  if(fragments %>% filter(type != "exon") %>% nrow() > 0){
    for(i in 1:nrow(fragments %>% filter(type != "exon"))){
      region_domain = fragments %>% 
        filter(type != "exon") %>%
        .[i,]
      region_domain_mod = region_domain %>%
        mutate(start = start * 3,
               end = end * 3)
      
      ### Add domain colored-filled-rectangles, start and stop aa 
      if (region_domain$type == "domain") {
        graph_gene = graph_gene +
          annotate("rect", xmin = region_domain_mod$start, xmax = region_domain_mod$end, ymin = 0, ymax = cds_seq_size/100,
                   fill = colors[1]) #+
      }
      
      
      if (region_domain$type == "zn_fingers") {
        graph_gene = graph_gene +
          annotate("rect", xmin = region_domain_mod$start, xmax = region_domain_mod$end, ymin = 0, ymax = cds_seq_size/100,
                   fill = colors[3]) #+
      }
      
      if (region_domain$type == "transmembrane") {
        graph_gene = graph_gene +
          annotate("rect", xmin = region_domain_mod$start, xmax = region_domain_mod$end, ymin = 0, ymax = cds_seq_size/100,
                   fill = colors[5]) #+
      }
      
      if (region_domain$type == "coiledCoil") {
        graph_gene = graph_gene +
          annotate("rect", xmin = region_domain_mod$start, xmax = region_domain_mod$end, ymin = 0, ymax = cds_seq_size/100,
                   fill = colors[2]) #+
      }
      
      if (region_domain$type == "NHL_YD_WD_repeats") {
        graph_gene = graph_gene +
          annotate("rect", xmin = region_domain_mod$start, xmax = region_domain_mod$end, ymin = 0, ymax = cds_seq_size/100,
                   fill = colors[4]) #+
      }
    }
  }
  
  
  mutation = gene_for_scheme %>%
    dplyr::filter(gene_name == g) %>%
    mutate(hjust = case_when(adjust == "leftleft" ~ 1.5,
                             adjust == "left" ~ 1,
                             adjust == "right" ~ 0,
                             TRUE ~ 0.5))
  
  if(any(!is.na(mutation$cDNA_mut_label))){
    for (m in 1:nrow(mutation)) {
      mut = mutation[m,]
      
      hjust = mut$hjust
      color = mut$color
      
      graph_gene = graph_gene +
        annotate("segment", x = mut$cDNA_position, xend = mut$cDNA_position, 
                 y = cds_seq_size/100*2.5, yend = cds_seq_size/100*4.5 + cds_seq_size/100*mut$correction, size = 1.3,
                 colour = color) +
        annotation_custom(grob=circleGrob(r=unit(3,"npc"), 
                                          gp = gpar(col = color, lty = 2, fill = color)), 
                          xmin = mut$cDNA_position - cds_seq_size/100*0.5, 
                          xmax = mut$cDNA_position + cds_seq_size/100*0.5, 
                          ymin=cds_seq_size/100*4.4 + cds_seq_size/100*mut$correction, 
                          ymax=cds_seq_size/100*4.52 + cds_seq_size/100*mut$correction) +
        annotate("text", x = mut$cDNA_position, y = cds_seq_size/100*6.4 + cds_seq_size/100*mut$correction, 
                 label = mut$cDNA_mut_label, hjust = hjust, vjust = 0.5, size = 5) +
        annotate("text", x = mut$cDNA_position, y = cds_seq_size/100*5.4 + + cds_seq_size/100*mut$correction, 
                 label = mut$protein_mut_label, hjust = hjust, vjust = 0.5, size = 5)
    }
    
    for(i in 1:nrow(fragments %>% filter(type == "exon"))){
      exon = fragments %>% 
        filter(type == "exon") %>%
        .[i,]
      
      ### Add exon and numbers
      graph_gene = graph_gene +
        annotate("segment", x = exon$start, xend = exon$start, 
                 y = cds_seq_size/100*1.5, yend = cds_seq_size/100*2.5,
                 colour = "black") +
        annotate("text", x = ((exon$end - exon$start + 1) / 2 + exon$start), y = cds_seq_size/100*2, label = exon$name,
                 size = 5, fontface = 2)
      
      
      ### Add exons start position (except first exon)
      
      graph_gene = graph_gene +
        annotate("text", x = exon$start, y = cds_seq_size/100*2.53, label = exon$start, hjust = 0, vjust = 1, size = 4, angle = 90)
    }
    
  }
  
  ggsave(paste0("results/scheme_genes/graph_", g, ".png"), width = 20, height = 2.7, dpi = 600, graph_gene) # for genes with variants in two levels
}  


## Selected gene TENM1

for (g in c("TENM1")){
  print(g)
  
  tbl = gene_for_scheme %>%
    filter(gene_name == g)
  
  transcript = tbl$RefSeq_transcript_id %>% 
    unique()
  
  protein_id = tbl$Uniprot_id %>%
    unique()
  
  exons = exons_coordinates %>%
    filter(external_gene_name == g) %>%
    mutate(name = as.character(name)) %>%
    dplyr::select(name, start, end, type)
  
  domains = domains_list %>%
    filter(protein == protein_id) %>%
    dplyr::select(name, start, end, type) 
  
  zn_fingers = zn_fingers_list %>%
    filter(protein == protein_id) %>%
    dplyr::select(name, start, end, type)
  
  repeats = repeats_list %>%
    filter(protein == protein_id) %>%
    dplyr::select(name, start, end, type)
  
  coiledCoil = coiledCoil_list %>%
    filter(protein == protein_id) %>%
    dplyr::select(name, start, end, type)
  
  tm = tm_list %>%
    filter(protein == protein_id) %>%
    dplyr::select(name, start, end, type)
  
  cds_seq_size = max(exons$end) - min(exons$start) + 1 - 3 # without the stop codon
  protein_size = cds_seq_size / 3
  
  fragments = rbind(exons, domains, zn_fingers, repeats, coiledCoil, tm) %>%
    mutate(position = end)
  
  position_max = max(exons_coordinates$end) + 200
  
  graph_gene = fragments %>%
    ggplot(aes(position)) +
    
    ### Add transcript and protein IDs
    labs(title = substitute(italic(x)/y/z, list(x = g, y = transcript, z = protein_id))) +
    coord_cartesian(clip = "off") +
    xlim(0, position_max) + ###
    theme_classic() +
    theme(
      plot.title = element_text(size = 20),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    
    ### Add transcript rectangle and first and last positions
    annotate("rect", xmin = 1, xmax = cds_seq_size, ymin = cds_seq_size/100*1.5, ymax = cds_seq_size/100*2.5,
             color = "black", fill = "white") +
    annotate("text", x = cds_seq_size + cds_seq_size*0.02, y = cds_seq_size/100*2, label = paste0(cds_seq_size, " bp"),
             hjust = 0) +
    
    ### Add protein rectangle and first and last positions
    annotate("rect", xmin = 1, xmax = cds_seq_size, ymin = 0, ymax = cds_seq_size/100,
             color = "black", fill = "white") +
    annotate("text", x = cds_seq_size + cds_seq_size*0.02, y = cds_seq_size/100/2, label = paste0(protein_size, " aa"),
             hjust = 0)
  
  
  if(fragments %>% filter(type != "exon") %>% nrow() > 0){
    for(i in 1:nrow(fragments %>% filter(type != "exon"))){
      region_domain = fragments %>% 
        filter(type != "exon") %>%
        .[i,]
      region_domain_mod = region_domain %>%
        mutate(start = start * 3,
               end = end * 3)
      

      ### Add domain colored-filled-rectangles, start and stop aa 
      if (region_domain$type == "domain") {
        graph_gene = graph_gene +
          annotate("rect", xmin = region_domain_mod$start, xmax = region_domain_mod$end, ymin = 0, ymax = cds_seq_size/100,
                   fill = colors[1]) #+
      }
      
      
      if (region_domain$type == "zn_fingers") {
        graph_gene = graph_gene +
          annotate("rect", xmin = region_domain_mod$start, xmax = region_domain_mod$end, ymin = 0, ymax = cds_seq_size/100,
                   fill = colors[3]) #+
      }
      
      if (region_domain$type == "transmembrane") {
        graph_gene = graph_gene +
          annotate("rect", xmin = region_domain_mod$start, xmax = region_domain_mod$end, ymin = 0, ymax = cds_seq_size/100,
                   fill = colors[5]) #+
      }
      
      if (region_domain$type == "coiledCoil") {
        graph_gene = graph_gene +
          annotate("rect", xmin = region_domain_mod$start, xmax = region_domain_mod$end, ymin = 0, ymax = cds_seq_size/100,
                   fill = colors[2]) #+
      }
      
      if (region_domain$type == "NHL_YD_WD_repeats") {
        graph_gene = graph_gene +
          annotate("rect", xmin = region_domain_mod$start, xmax = region_domain_mod$end, ymin = 0, ymax = cds_seq_size/100,
                   fill = colors[4]) #+
      }
    }
  }
  
  
  mutation = gene_for_scheme %>%
    dplyr::filter(gene_name == g) %>%
    mutate(hjust = case_when(adjust == "leftleft" ~ 1.5,
                             adjust == "left" ~ 1,
                             adjust == "right" ~ 0,
                             TRUE ~ 0.5))
  
  if(any(!is.na(mutation$cDNA_mut_label))){
    for (m in 1:nrow(mutation)) {
      mut = mutation[m,]
      
      hjust = mut$hjust
      color = mut$color
      
      graph_gene = graph_gene +
        annotate("segment", x = mut$cDNA_position, xend = mut$cDNA_position, 
                 y = cds_seq_size/100*2.5, yend = cds_seq_size/100*4.5 + cds_seq_size/100*mut$correction, size = 1.3,
                 colour = color) +
        annotation_custom(grob=circleGrob(r=unit(3,"npc"), 
                                          gp = gpar(col = color, lty = 2, fill = color)), 
                          xmin = mut$cDNA_position - cds_seq_size/100*0.5, 
                          xmax = mut$cDNA_position + cds_seq_size/100*0.5, 
                          ymin=cds_seq_size/100*4.4 + cds_seq_size/100*mut$correction, 
                          ymax=cds_seq_size/100*4.52 + cds_seq_size/100*mut$correction) +
        annotate("text", x = mut$cDNA_position, y = cds_seq_size/100*6.4 + cds_seq_size/100*mut$correction, 
                 label = mut$cDNA_mut_label, hjust = hjust, vjust = 0.5, size = 5) +
        annotate("text", x = mut$cDNA_position, y = cds_seq_size/100*5.4 + + cds_seq_size/100*mut$correction, 
                 label = mut$protein_mut_label, hjust = hjust, vjust = 0.5, size = 5)
    }
    
    for(i in 1:nrow(fragments %>% filter(type == "exon"))){
      exon = fragments %>% 
        filter(type == "exon") %>%
        .[i,]
      
      ### Add exon and numbers
      graph_gene = graph_gene +
        annotate("segment", x = exon$start, xend = exon$start, 
                 y = cds_seq_size/100*1.5, yend = cds_seq_size/100*2.5,
                 colour = "black") +
        annotate("text", x = ((exon$end - exon$start + 1) / 2 + exon$start), y = cds_seq_size/100*2, label = exon$name,
                 size = 5, fontface = 2)
      
      
      ### Add exons start position (except first exon)
      
      graph_gene = graph_gene +
        annotate("text", x = exon$start, y = cds_seq_size/100*2.53, label = exon$start, hjust = 0, vjust = 1, size = 4, angle = 90)
    }
    
  }
  
  ggsave(paste0("results/scheme_genes/graph_", g, ".png"), width = 20, height = 2.7, dpi = 600, graph_gene)
  
}  


## Selected genes validation

for (g in unique(gene_for_scheme_val$gene_name)){
  
  tbl = gene_for_scheme_val %>%
    filter(gene_name == g)
  
  transcript = tbl$RefSeq_transcript_id %>% 
    unique()
  
  protein_id = tbl$Uniprot_id %>%
    unique()
  
  exons = exons_coordinates_val %>%
    filter(external_gene_name == g) %>%
    mutate(name = as.character(name)) %>%
    dplyr::select(name, start, end, type)
  
  domains = domains_list_val %>%
    filter(protein == protein_id) %>%
    dplyr::select(name, start, end, type)
  
  repeats = repeats_list_val %>%
    filter(protein == protein_id) %>%
    dplyr::select(name, start, end, type)
  
  tm = tm_list_val %>%
    filter(protein == protein_id) %>%
    dplyr::select(name, start, end, type)
  
  cds_seq_size = max(exons$end) - min(exons$start) + 1 - 3 # without the stop codon
  protein_size = cds_seq_size / 3
  
  fragments = rbind(exons, domains, repeats, tm) %>%
    mutate(position = end)
  
  graph_gene = fragments %>%
    ggplot(aes(position)) +
    
    ### Add transcript and protein IDs
    labs(title = substitute(italic(x)/y/z, list(x = g, y = transcript, z = protein_id))) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 20),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    
    ### Add transcript rectangle and first and last positions
    annotate("rect", xmin = 1, xmax = cds_seq_size, ymin = cds_seq_size/100*1.5, ymax = cds_seq_size/100*2.5,
             color = "black", fill = "white") +
    annotate("text", x = cds_seq_size + cds_seq_size*0.02, y = cds_seq_size/100*2, label = paste0(cds_seq_size, " bp"),
             hjust = 0) +
    
    ### Add protein rectangle and first and last positions
    annotate("rect", xmin = 1, xmax = cds_seq_size, ymin = 0, ymax = cds_seq_size/100,
             color = "black", fill = "white") +
    annotate("text", x = cds_seq_size + cds_seq_size*0.02, y = cds_seq_size/100/2, label = paste0(protein_size, " aa"),
             hjust = 0)
  
  
  if(fragments %>% filter(type != "exon") %>% nrow() > 0){
    for(i in 1:nrow(fragments %>% filter(type != "exon"))){
      region_domain = fragments %>% 
        filter(type != "exon") %>%
        .[i,]
      region_domain_mod = region_domain %>%
        mutate(start = start * 3,
               end = end * 3)
      

      ### Add domain colored-filled-rectangles, start and stop aa 
      
      if (region_domain$type == "domain") {
        graph_gene = graph_gene +
          annotate("rect", xmin = region_domain_mod$start, xmax = region_domain_mod$end, ymin = 0, ymax = cds_seq_size/100,
                   fill = colors[1])# +
      }
      
      if (region_domain$type == "transmembrane") {
        graph_gene = graph_gene +
          annotate("rect", xmin = region_domain_mod$start, xmax = region_domain_mod$end, ymin = 0, ymax = cds_seq_size/100,
                   fill = colors[5])# +
      }
      
      if (region_domain$type == "ANK_repeats") {
        graph_gene = graph_gene +
          annotate("rect", xmin = region_domain_mod$start, xmax = region_domain_mod$end, ymin = 0, ymax = cds_seq_size/100,
                   fill = colors[4]) #+
      }
    }
  }
  
  
  mutation = gene_for_scheme_val %>%
    dplyr::filter(gene_name == g) %>%
    mutate(hjust = case_when(adjust == "leftleft" ~ 1.5,
                             adjust == "left" ~ 1,
                             adjust == "right" ~ 0,
                             TRUE ~ 0.5))
  
  if(any(!is.na(mutation$cDNA_mut_label))){
    for (m in 1:nrow(mutation)) {
      mut = mutation[m,]
      
      hjust = mut$hjust
      color = mut$color
      
      
      graph_gene = graph_gene +
        annotate("segment", x = mut$cDNA_position, xend = mut$cDNA_position, 
                 y = cds_seq_size/100*2.5, yend = cds_seq_size/100*4.5 + cds_seq_size/100*mut$correction, size = 1.3,
                 colour = color) +
        annotation_custom(grob=circleGrob(r=unit(3,"npc"), 
                                          gp = gpar(col = color, lty = 2, fill = color)), 
                          xmin = mut$cDNA_position - cds_seq_size/100*0.5, 
                          xmax = mut$cDNA_position + cds_seq_size/100*0.5, 
                          ymin=cds_seq_size/100*4.4 + cds_seq_size/100*mut$correction, 
                          ymax=cds_seq_size/100*4.52 + cds_seq_size/100*mut$correction) +
        annotate("text", x = mut$cDNA_position, y = cds_seq_size/100*6.4 + cds_seq_size/100*mut$correction, 
                 label = mut$cDNA_mut_label, hjust = hjust, vjust = 0.5, size = 5) +
        annotate("text", x = mut$cDNA_position, y = cds_seq_size/100*5.4 + + cds_seq_size/100*mut$correction, 
                 label = mut$protein_mut_label, hjust = hjust, vjust = 0.5, size = 5)
    }
    
    for(i in 1:nrow(fragments %>% filter(type == "exon"))){
      exon = fragments %>% 
        filter(type == "exon") %>%
        .[i,]
      
      ### Add exon and numbers
      graph_gene = graph_gene +
        annotate("segment", x = exon$start, xend = exon$start, 
                 y = cds_seq_size/100*1.5, yend = cds_seq_size/100*2.5,
                 colour = "black") +
        annotate("text", x = ((exon$end - exon$start + 1) / 2 + exon$start), y = cds_seq_size/100*2, label = exon$name,
                 size = 5, fontface = 2)
      
      
      ### Add exons start position (except first exon)
      
      graph_gene = graph_gene +
        annotate("text", x = exon$start, y = cds_seq_size/100*2.53, label = exon$start, hjust = 0, vjust = 1, size = 4, angle = 90)
    }
    
  }
  
  ggsave(paste0("results/scheme_genes/graph_", g, ".png"), width = 18, height = 2.4, dpi = 600, graph_gene)

}  
