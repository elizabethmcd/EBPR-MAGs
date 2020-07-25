library(tidyverse)

# Runella annotations

runella_kofam <- read.delim("results/2013_binning/annotations/runella/3300009517-bin.35-annotations.txt", sep="\t", header=FALSE)
colnames(runella_kofam) <- c("locus_tag", "KO", "Kofam_annotation")
runella_prokka <- read.delim("results/2013_binning/annotations/runella/3300009517-bin.35-prokka-annotations.tsv", sep="\t", header=FALSE)
colnames(runella_prokka) <- c("locus_tag", "prokka_annotation")
merged_annotations <- left_join(runella_prokka, runella_kofam)

# Runella counts

runella_counts <- read.delim("results/2013_binning/annotations/runella/3300009517-bin.35-normalized-counts.tsv", sep=";", header=FALSE) %>% select(-V8, -V9)
colnames(runella_counts) <- c("locus_tag", "an1", "an2", "an3", "aer1", "aer2", "aer3")

# Merged counts and annotations
runella_table <- left_join(runella_counts, merged_annotations)
write.csv(runella_table, "results/2013_binning/annotations/runella/runella-counts-table.csv", quote=FALSE, row.names = FALSE)
