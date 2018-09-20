#! /usr/local/bin/r

# Packages
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

options( warn = -1 )

# Read in datasets
all_ani = read.delim(args[1], sep="\t", header=FALSE)
checkm_results = read.delim(args[2], sep=' ', header=FALSE)
colnames(all_ani) = c("bin", "bin2", "ANI1", "ANI2", "AF1", "AF2")
colnames(checkm_results) = c("bin", "classification", "size", "completeness", "redundancy")

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
colnames(full) = c("bin1","bin2","ANI1","ANI2","AF1","AF2","classf1","size1","completeness1","redundancy1","classf2","size2","completeness2","redundancy2")

# Subset the full dataset to only greater than 90% complete to get rid of spurious small things
greater90both = full %>% filter(completeness1 > 90 & completeness2 > 90)
greater90comp = full %>% filter(completeness1 > 90 | completeness2 > 90)

# Parse out the identical ones
identical = greater90comp %>% filter(ANI1 > 99 | ANI2 > 99)

# write out file
write.table(identical, file="identical-bins.txt", sep="\t")
