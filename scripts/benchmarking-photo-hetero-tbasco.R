library(tidyverse)
library(tximport)
library(stringr)

# count files 
dir <- "kallisto_files"
samples <- read.table(file.path(dir, "sra-list.txt"), header=FALSE)
files <- file.path(dir, samples$V1, "abundance.h5")
names(files) <- paste0("sample", 1:20)
txi.kallisto <- tximport(files, type="kallisto", txOut=TRUE)
counts <- as.data.frame(txi.kallisto)
finalcounts <- rownames_to_column(counts, var="locus_tag")

# KEGG annotations 
ko <- read.delim("genomes/ko_annotations_table.tsv", sep="\t", header=FALSE)
colnames(ko) <- c("genome", "locus_tag", "annotation")

# merge
ko_counts <- left_join(ko, finalcounts)

# tables of raw vs TPM normalized tables
raw <- ko_counts[,(1:23)]
norm <- ko_counts[,c(1:3,24:43)]

# TbasCO formatted input tables
# order of locus_tag, sample counts, annotation, bin
tbasco_raw <- raw[,c(2,4:23,3,1)]
tbasco_norm <- norm[,c(2,4:23,3,1)]

# total raw counts across all samples for each genome
test <- tbasco_raw
sumCounts <- aggregate(test[2:21], list(test$genome), sum)
sumCounts$total_counts <- rowSums(sumCounts[,(2:21)])

# write out tables
write_delim(tbasco_raw, "tables/raw_benchmarking_counts_annots.tsv", delim=";")

# normalized counts within McClure et al. for comparison
mcclure = read.csv("tables/mcclure_normalized.tsv")
new = mcclure %>% separate(GeneID, into = c("genome", "ID"), sep=3)
sumMcclure <- aggregate(new[3:27], list(new$genome), sum)
sumMcclure$total_counts <- rowSums(sumMcclure[,(2:26)])
