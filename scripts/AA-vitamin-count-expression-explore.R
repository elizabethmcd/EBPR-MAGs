library(tidyverse)

# metadata
metadata <- read.csv("results/R1R2-EBPR-MAGs-table.csv")

# vitamin annotations
vitamin_list <- read.csv("results/2013_binning/annotations/vitamin-list.csv", header = FALSE)
colnames(vitamin_list) <- c("ko", "pathway")
vitamin_annotations <- read.delim("results/2013_binning/annotations/R1R2-vitamin-annotations.txt", sep = "\t", header=FALSE)
colnames(vitamin_annotations) <- c("locus_tag", "ko", "annotation")
vitamin_table <- left_join(vitamin_annotations, vitamin_list)
vitamin_table$Bin <- gsub("\\_.*", "", vitamin_table$locus_tag)
vitamin_results <- vitamin_table %>% select(Bin, ko, annotation, pathway) %>% arrange(Bin)
vitamin_info <- left_join(vitamin_results, metadata) %>% select(Bin, Code, ko, annotation, pathway)

# amino acid annotations
amino_list <- read.csv("results/2013_binning/annotations/amino_acid_list.csv", header=FALSE)
colnames(amino_list) <- c("ko", "pathway")
amino_annotations <- read.delim("results/2013_binning/annotations/R1R2-amino-acid-annotations.txt", sep="\t", header=FALSE)
colnames(amino_annotations) <- c("locus_tag", "ko", "annotation")
amino_table <- left_join(amino_annotations, amino_list)
amino_table$Bin <- gsub("\\_.*", "", amino_table$locus_tag)
amino_results <- amino_table %>% select(Bin, ko, annotation, pathway) %>% arrange(Bin)
amino_info <- left_join(amino_results, metadata) %>% select(Bin, Code, ko, annotation, pathway)

# IIA tables
IIA_amino_table <- amino_info %>% filter(Code=="CAPIIA") %>% arrange(pathway)
IIA_vitamin_table <- vitamin_info %>% filter(Code=="CAPIIA") %>% arrange(pathway)

# IA tables
IA_amino_table <- amino_info %>% filter(Code=="CAPIA") %>% arrange(pathway)
IA_vitamin_table <- vitamin_info %>% filter(Code=="CAPIA") %>% arrange(pathway)

# Totals
vitamin_totals <- vitamin_info %>% group_by(Code) %>% count() %>% arrange(n)
amino_totals <- amino_info %>% group_by(Code) %>% count()

# vitamins
biotin_count <- vitamin_info %>% group_by(Code) %>% filter(pathway=='biotin ') %>% count() %>% arrange(n)
cobalamin_count <- vitamin_info %>% group_by(Code) %>% filter(pathway=="cobalamin") %>% count() %>% arrange(n)
THF_count <- vitamin_info %>% group_by(Code) %>% filter(pathway=="THF") %>% count() %>% arrange(n)
thiamine_count <- vitamin_info %>% group_by(Code) %>% filter(pathway=="thiamine") %>% count() %>% arrange(n)

##################################################
# Expression data added
##################################################
vitamin_expression <- read.delim("results/2013_binning/annotations/vitamin-expression-table.csv", sep=";", header=FALSE)
amino_expression <- read.delim("results/2013_binning/annotations/amino-expression-table.csv", sep=';', header=FALSE)
names <- c("locus_tag", "t1", "t2", "t3", "t4", "t5", "t6", "ko", "Bin")
colnames(vitamin_expression) <- c(names)
colnames(amino_expression) <- c(names)

# vitamin table
vitamin_pathways <- left_join(vitamin_expression, vitamin_list)
vitamin_expression_table <- left_join(vitamin_pathways, metadata) %>% 
  select(Bin, Code, t1, t2, t3, t4, t5, t6, ko, pathway) %>% 
  arrange(Code)

CAP_vitmain_expression <- vitamin_expression_table %>% filter(Code=="CAPIIA"| Code =='CAPIA') %>% arrange(Code,pathway)

# amino acid table
amino_pathways <- left_join(amino_expression, amino_list)
amino_expression_table <- left_join(amino_pathways, metadata) %>% 
  select(Bin, Code, t1, t2, t3, t4, t5, t6, ko, pathway) %>% 
  arrange(Code)

CAP_amino_expression <- amino_expression_table %>% filter(Code=="CAPIIA" | Code=="CAPIA") %>% arrange(Code,pathway)

# top 20 expressed genomes 

top_20 <- c("CAPIIA",
            "RUN1",
            "BAC3",
            "CAULO1",
            "HYPHO1",
            "PSEUDO1",
            "RHODO1",
            "CHIT1",
            "FLAVO1",
            "GEMMA1",
            "TET2",
            "CAPIA",
            "RAM1",
            "OBS1",
            "TET1",
            "ALPHA1",
            "RUBRI1",
            "LEAD1",
            "ZOO1")

top_15 <- c("CAPIIA",
            "RUN1",
            "BAC3",
            "CAULO1",
            "HYPHO1",
            "PSEUDO1",
            "RHODO1",
            "CHIT1",
            "FLAVO1",
            "GEMMA1",
            "TET2",
            "CAPIA",
            "RAM1",
            "OBS1",
            "TET1")

top_20_vitamin_expression <- vitamin_expression_table %>% 
  filter(Code %in% top_20) %>% 
  arrange(Code,pathway)

top_20_amino_expression <- amino_expression_table %>% 
  filter(Code %in% top_20) %>% 
  arrange(Code, pathway)

top_15_vitamin_expression <- vitamin_expression_table %>% 
  filter(Code %in% top_15) %>% 
  arrange(Code, pathway)

top_15_amino_expression <- amino_expression_table %>% 
  filter(Code %in% top_15) %>% 
  arrange(Code, pathway)

# Totals for top 15
vitamin_top_15_totals <- vitamin_info %>% 
  filter(Code %in% top_15) %>% 
  group_by(Code) %>% 
  count() %>% 
  arrange(n)

amino_top_15_totals <- amino_info %>% 
  filter(Code %in% top_15) %>% 
  group_by(Code) %>% 
  count() %>% 
  arrange(n)

amino_top_15_descriptions <- amino_info %>% 
  filter(Code %in% top_15) %>% 
  group_by(Code) %>% 
  arrange(Code, pathway)

amino_top_15_totals

write.csv(amino_top_15_descriptions, "results/2013_binning/annotations/manual_annotations/top15_amino_descriptions.csv", quote=FALSE, row.names = FALSE)
