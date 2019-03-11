library(tximport)
library(tximportData)
library(readr)
library(tibble)
library(dplyr)

# kallisto files
dir <- "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/kallisto_results"
samples <- read.table(file.path(dir, "samples.txt"), header=TRUE)
files <- file.path(dir, samples$experiment, "abundance.h5")
names(files) <- paste0("sample", 1:6)
txi.kallisto <- tximport(files, type="kallisto", txOut = TRUE)
counts <- as.data.frame(txi.kallisto)
finalcounts <- rownames_to_column(counts, var="ID")
counttable <- finalcounts[, c(1,8:13)]
write_delim(counttable, "~/Desktop/kallisto_test.tsv", delim="\t")

# merge with KO annotations and metadata on genomes
ko <- read.delim("~/Desktop/ebpr-ko-annotations.txt", sep="\t", header=FALSE)
metadata <- read.delim("/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/results/refined-bins-derep-checkm-stats.txt", header=FALSE, sep=" ")
orfs <- read.delim("~/Desktop/ebpr-orfs.txt", header=FALSE, sep="\t")
colnames(orfs) <- c("Genome", "orfs")
colnames(ko) <- c("Genome", "ID", "KO")
ko_counts <- left_join(counttable, ko)
sumCounts <- aggregate(ko_counts[2:7], list(ko_counts$Genome), sum)
sumCounts$Mean <- rowMeans(sumCounts[,-1])
colnames(sumCounts)[1] <- "Genome"
countSizes <- left_join(sumCounts, orfs)
countSizes$orfs80 <- countSizes$orfs * .80
countSizes$cutoff_score <- countSizes$Mean / countSizes$orfs80

# merge with Prokka annotations, filtering out tRNA's for normalization purposes
prokka <- read.delim("~/Desktop/ebpr-prokka-annotations.txt", sep="\t", header=FALSE)
colnames(prokka) <- c("ID", "prokka_annotation", "size_bp", "accession")
counts_all_annots <- left_join(prokka, ko_counts)
ebpr_counts <- counts_all_annots[,c(1, 5:10, 12, 2:4, 11)]

# determined cutoffs, create countTable with those above cutoff for downstream analysis
above_cutoff <- c("3300009517-bin.1", "3300009517-bin.12", "3300009517-bin.13", "3300009517-bin.29", "3300009517-bin.3", "3300009517-bin.30", "3300009517-bin.31", "3300009517-bin.42", "3300009517-bin.47", "3300009517-bin.52", "3300009517-bin.6", "3300009517-bin.7", "3300026282-bin.4", "3300026283-bin.21", "3300026283-bin.28", "3300026284-bin.6", "3300026284-bin.9", "3300026287-bin.29", "3300026287-bin.4", "3300026288-bin.19", "3300026288-bin.32", "3300026288-bin.43", "3300026289-bin.24", "3300026289-bin.41", "3300026299-bin.22", "3300026299-bin.26", "3300026302-bin.10", "3300026302-bin.20", "3300026302-bin.25", "3300026302-bin.31", "3300026302-bin.32", "3300026302-bin.46", "3300026302-bin.47", "3300026302-bin.62", "3300026303-bin.42", "3300026303-bin.46")

ebpr_above_cutoff <- ebpr_counts %>% filter(Genome %in% above_cutoff)
ebpr_above_cutoff$size_kbp <- ebpr_above_cutoff$size_bp / 1000
ebpr_above_cutoff$rpk1 <- ebpr_above_cutoff$counts.sample1 / ebpr_above_cutoff$size_kbp
ebpr_above_cutoff$rpk2 <- ebpr_above_cutoff$counts.sample2 / ebpr_above_cutoff$size_kbp
ebpr_above_cutoff$rpk3 <- ebpr_above_cutoff$counts.sample3 / ebpr_above_cutoff$size_kbp
ebpr_above_cutoff$rpk4 <- ebpr_above_cutoff$counts.sample4 / ebpr_above_cutoff$size_kbp
ebpr_above_cutoff$rpk5 <- ebpr_above_cutoff$counts.sample5 / ebpr_above_cutoff$size_kbp
ebpr_above_cutoff$rpk6 <- ebpr_above_cutoff$counts.sample6 / ebpr_above_cutoff$size_kbp

PM1 <- sum(ebpr_above_cutoff$rpk1) / 1000000
PM2 <- sum(ebpr_above_cutoff$rpk2) / 1000000
PM3 <- sum(ebpr_above_cutoff$rpk3) / 1000000
PM4 <- sum(ebpr_above_cutoff$rpk4) / 1000000
PM5 <- sum(ebpr_above_cutoff$rpk5) / 1000000
PM6 <- sum(ebpr_above_cutoff$rpk6) / 1000000

ebpr_above_cutoff$TPM1 <- ebpr_above_cutoff$rpk1 / PM1
ebpr_above_cutoff$TPM2 <- ebpr_above_cutoff$rpk2 / PM2
ebpr_above_cutoff$TPM3 <- ebpr_above_cutoff$rpk3 / PM3
ebpr_above_cutoff$TPM4 <- ebpr_above_cutoff$rpk4 / PM4
ebpr_above_cutoff$TPM5 <- ebpr_above_cutoff$rpk5 / PM5
ebpr_above_cutoff$TPM6 <- ebpr_above_cutoff$rpk6 / PM6

tbasco_table <- ebpr_above_cutoff[,c(1,20:25,8,12)]
colnames(tbasco_table) <- c("Locus_Tag", "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6", "Annotation", "Bin")
raw_table <- ebpr_above_cutoff[,c(1:7,8,12)]
colnames(raw_table) <- c("Locus_Tag", "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6", "Annotation", "Bin")

# calculatng TPM for all counts for normalization by sample (by genome size?)
# divide read counts by size of gene in kb = reads per kilobase (RPK)
# add up all RPK values in a sample and divide by 1,000,000 = per million scaling factor
# divide RPK values by per million scaling factor = TPM

write_delim(countSizes, "~/Desktop/ebpr-transcriptomes-cutoffs.tsv", delim="\t")
write_delim(tbasco_table, "~/Desktop/tbasco-sample-data.csv", delim=";")
write_delim(raw_table, "~/Desktop/raw-ebpr-counts-data.csv", delim=";")
