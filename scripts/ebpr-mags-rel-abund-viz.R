# Analyzing EBPR Refined Bins Relative Abundance
library(tidyverse)
library(reshape2)

# Files 
covg = read.delim("results/EBPR-bins-final-mapped-clean.coverage.txt", sep="\t", header=FALSE)
taxonomy = read.delim("results/gtdbtk.bac120.classification_pplacer.txt", sep="\t", header=FALSE)

# cleanup coverage file to get relabund table
covgprep = covg %>% select(c("V2", "V3", "V8"))
covtable = as.data.frame(spread(covgprep, key="V2", value="V8", fill=0))
rownames(covtable) = covtable[,1]
covtable = subset(covtable, select = -c(V3))
covtab = as.data.frame(lapply(covtable, as.numeric))
binnames = colnames(covtable)
samplenames = rownames(covtable)
colnames(covtab) = binnames
rownames(covtab) = samplenames
rowSums(covtab) / 57
reltab = covtab; for(i in seq_len(nrow(reltab))) reltab[i, ] = reltab[i, ] / rowSums(covtab)

# cleanup taxonomy file 
colnames(taxonomy) = c("Bin", "Taxonomy")
tax.clean = separate(taxonomy, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern="d__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern="p__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern="c__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern="o__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern="f__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern="g__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern="s__", replacement=""))
row.names(tax.clean) = tax.clean[,1]
tax.clean = tax.clean[,-1]

# Sum abundance across samples for a rank abundance for all genomes
taxtot = rownames_to_column(tax.clean, "Bin")
totalrel = as.data.frame(colSums(reltab) / 10)
totalcov = as.data.frame(colSums(covtab) / 10)
totalrel = rownames_to_column(totalrel, "Bin")
totalcov = rownames_to_column(totalcov, "Bin")
colnames(totalrel) = c("Bin", "Abundance")
colnames(totalcov) =c("Bin", "Abundance")
ebpr = left_join(totalrel, taxtot)
ebpr_cov = left_join(totalcov, taxtot)
ebpr_ordered = ebpr[order(ebpr$Abundance, decreasing=T), ]
ebpr_cov_ord = ebpr_cov[order(ebpr_cov$Abundance, decreasing=T), ]

# write out joined ebpr dataframe
write.csv(file ="results/ebpr-rank-abundance-classifications.csv", ebpr)
write.csv(file = "results/ebpr-rank-cov-classifications.csv", ebpr_cov_ord)

# read back in with manually changed classifications
ebpr_maned = read.csv("results/ebpr-rank-abundance-classifications.csv", header=TRUE)
ebpr_coved = read.csv("results/ebpr-rank-cov-classifications.csv", header=TRUE)

# plot rank abundance
relabundplot <- ggplot(ebpr_maned, aes(x=reorder(Bin,-Abundance), y=Abundance, fill=Highest_Classf)) + geom_col() + labs(x="Genome", y="Average Abundance Across all 10 Samples") + theme_classic() + scale_fill_brewer(palette = "Paired") 

covplot <- ggplot(ebpr_coved, aes(x=reorder(Bin,-Abundance), y=Abundance, fill=Highest_Classf)) + geom_col() + labs(x="Genome", y="Average Abundance Across all 10 Samples") + theme_classic() + scale_fill_brewer(palette = "Paired")

phylaplot <- ggplot(ebpr_maned, aes(x=reorder(Bin,-Abundance), y=Abundance, fill=Phylum)) + geom_col() + labs(x="Genome", y="Average Abundance Across all 10 Samples") + theme_classic() + scale_fill_brewer(palette = "Paired") 

# bubble plot of time series points
relabund_prep = rownames_to_column(reltab, "meta")
relhigh = gather(relabund_prep, key="meta", value="Bin")
colnames(relhigh) = c("Bin", "Abundance")
relhigh$meta = rep_len(samplenames, length.out=57)
relhighnames = left_join(relhigh, taxonomy)
relhigh.clean = separate(relhighnames, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
relhigh.clean = as.data.frame(sapply(relhigh.clean, gsub, pattern="d__", replacement=""))
relhigh.clean = as.data.frame(sapply(relhigh.clean, gsub, pattern="p__", replacement=""))
relhigh.clean = as.data.frame(sapply(relhigh.clean, gsub, pattern="c__", replacement=""))
relhigh.clean = as.data.frame(sapply(relhigh.clean, gsub, pattern="o__", replacement=""))
relhigh.clean = as.data.frame(sapply(relhigh.clean, gsub, pattern="f__", replacement=""))
relhigh.clean = as.data.frame(sapply(relhigh.clean, gsub, pattern="g__", replacement=""))
relhigh.clean = as.data.frame(sapply(relhigh.clean, gsub, pattern="s__", replacement=""))

ebpr_manames = ebpr_maned %>% select(c("Bin", "Highest_Classf"))
reltimetable = left_join(relhigh.clean, ebpr_manames)

ggsave(filename="ebpr-rank-abundance-plot.png", plot=relabundplot, width=40, height=20, units=c("cm"))
