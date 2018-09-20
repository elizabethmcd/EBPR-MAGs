#! /usr/local/bin/r

# Packages
library(dplyr)
library(data.table)
library(igraph)

options( warn = -1 ) 
# Read in datasets
all_ani = read.delim("results/EBPR-BINS.all.ani.out.cleaned", sep="\t", header=FALSE)
checkm_results = read.delim("results/EBPR-bins-checkm-results.txt", sep=' ', header=FALSE)
genome_stats = read.delim("results/all-bins-stats.txt", sep="\t", header=FALSE)

# Cleanup genome stats file
genome_stats_clean = genome_stats %>% filter(V1 != "n_scaffolds")
genome_stats_red = genome_stats_clean %>% select(c(V9,V20))
genome_stats_red$V20 = gsub("/Volumes/mcmahonlab/home/emcdaniel/EBPR-projects/EBPR-Flanking-Bins/all-bins-coded/", "", genome_stats_red$V20)
genome_stats_red$V20 = gsub(".fna","",genome_stats_red$V20)

# Column names
colnames(all_ani) = c("bin", "bin2", "ANI1", "ANI2", "AF1", "AF2")
colnames(checkm_results) = c("bin", "classification", "size", "completeness", "redundancy")
colnames(genome_stats_red) = c("contigL50", "bin")
checkm_results$size = checkm_results$size/1000000

# get rid of file extensions in all_ani file
all_ani$bin = gsub(".fna-genes.fna", "", all_ani$bin)
all_ani$bin2 = gsub(".fna-genes.fna", "", all_ani$bin2)

# Subsetting datasets for looking at alignment fraction
nonzeros = all_ani %>% filter(ANI1 != 0 & ANI2 != 0 & AF1 !=0 & AF2 !=0)
greater3 = nonzeros %>% filter(AF1 > 0.30 & AF2 > 0.30)
greater5 = nonzeros %>% filter(AF1 > 0.50 & AF2 > 0.50)

# Merge datasets and change column names 
new = left_join(greater5, checkm_results, by.y="bin")
colnames(new) = c("bin1", "bin", "ANI1", "ANI2", "AF1", "AF2", "classf1", "size1", "completeness1", "redundancy1")
full = left_join(new, checkm_results, by.y="bin")
colnames(full) = c("bin","bin2","ANI1","ANI2","AF1","AF2","classf1","size1","completeness1","redundancy1","classf2","size2","completeness2","redundancy2")
stats1 = left_join(full, genome_stats_red, by.y="bin")
colnames(stats1) = c("bin1","bin","ANI1","ANI2","AF1","AF2","classf1","size1","completeness1","redundancy1","classf2","size2","completeness2","redundancy2", "contigL50.1")
completedf = left_join(stats1, genome_stats_red, by.y="bin")
colnames(completedf) = c("bin1","bin2","ANI1","ANI2","AF1","AF2","classf1","size1","completeness1","redundancy1","classf2","size2","completeness2","redundancy2", "contigL50.1", "contigL50.2")

# Subset the full dataset to only greater than 90% complete to get rid of spurious small things
greater90both = completedf %>% filter(completeness1 > 90 & completeness2 > 90)
greater90comp = completedf %>% filter(completeness1 > 90 | completeness2 > 90)

# Parse out the identical ones
identical = greater90comp %>% filter(ANI1 > 99 | ANI2 > 99)

# Parse between 95 and 97%, within the same species, but also need to look at within the same sample 
actino = greater90comp %>% filter(classf1 == "'o__Actinomycetales'" & classf2 == "'o__Actinomycetales'")
xantho = greater90comp %>% filter(classf1 == "'f__Xanthomonadaceae'" & classf2 == "'f__Xanthomonadaceae'")
ctyo = greater90comp %>% filter(classf1 == "'o__Cytophagales'" & classf2 == "'o__Cytophagales'")

# Generate a graph from the dataframe to get high matching sets - idea from Ben
adjacency.graph = graph_from_data_frame(identical)
HMS.ID = data.frame(paste("HMS.", clusters(adjacency.graph)$membership, sep=""), names(clusters(adjacency.graph)$membership), stringsAsFactors=FALSE)
names(HMS.ID) = c("HMS", "bin")
HMS.ID = HMS.ID %>% arrange(HMS)

# write out files
write.table(identical, file="identical-bins.txt", sep="\t", row.names=FALSE)
write.csv(file="output.HMS", HMS.ID, row.names=FALSE)
