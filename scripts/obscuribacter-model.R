library(tidyverse)

# Obscuribacter annotations

obs_annotations <- read.delim("results/2013_binning/annotations/obscuribacter-kofam-sig-annotations.txt", sep="\t", header=FALSE)
colnames(obs_annotations) <- c("locus_tag", "KO", "annotation")
# Obscuribacter counts
counts <- read.delim("results/2013_transcriptomes/tables/obscuribacter_counts.tsv", sep=";", header=FALSE)
colnames(counts) <- c("locus_tag", "an1", "an2", "an3", "aer1", "aer2", "aer3", "KO", "genome")
obs_counts <- counts %>% select(-genome)

obs_prokka <- read.delim("results/2013_binning/annotations/obscuribacter_prokka_annotations.tsv", sep="\t", header=FALSE)
colnames(obs_prokka) <- c("locus_tag", "prokka_annotation")

# annotations and counts
obs_ko <- left_join(counts, obs_annotations)
obs_table <- left_join(obs_ko, obs_prokka)

write.csv(obs_table, "results/2013_binning/annotations/obscuri/obscuribacter_annotations_normalized_counts.csv", quote=FALSE, row.names = FALSE)
