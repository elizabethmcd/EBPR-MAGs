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
ebpr_names = ebpr_maned %>% select("Bin", "Highest_Classf")
relhighnames = left_join(relhigh, ebpr_names)
write.csv(file="results/relative-abundance-time-series-table.csv", relhighnames)

# relative abundance timeseries table 
reltimetable = read.csv("results/relative-abundance-time-series-table.csv", header=TRUE)
reltimetable$Abundance = as.numeric(as.character(reltimetable$Abundance))
reltimetable$Bin = gsub("-", "", reltimetable$Bin)
reltimetable = reltimetable %>% select(c("Bin", "Abundance", "Highest_Classf", "Sample"))

# only accumulibacter
accum = reltimetable %>% filter(Bin == "3300026282bin.4" | Bin == "3300026284bin.9")

accum_rel = ggplot(accum, aes(x=Sample, y=Bin)) + geom_point(aes(colour=Highest_Classf, size=Abundance)) + theme(panel.background = element_rect(fill="white"), legend.key=element_blank(), panel.grid.major=element_line(size=.5,colour="grey95", linetype="dotted"), panel.grid.minor=element_line(size=.5, colour="grey95", linetype="dotted"), axis.text.x=element_text(angle=60,size=10,vjust=.5), axis.text.y=element_text(size=7), axis.line=element_blank())

fullset = ggplot(reltimetable, aes(x=Sample, y=Bin)) + geom_point(aes(colour=Highest_Classf, size=Abundance)) + theme(panel.background = element_rect(fill="white"), legend.key=element_blank(), panel.grid.major=element_line(size=.5,colour="grey95", linetype="dotted"), panel.grid.minor=element_line(size=.5, colour="grey95", linetype="dotted"), axis.text.x=element_text(angle=60,size=10,vjust=.5), axis.text.y=element_text(size=7), axis.line=element_blank())

ggsave(filename="ebpr-rank-abundance-plot.png", plot=relabundplot, width=40, height=20, units=c("cm"))
ggsave(filename="figs/accum-rel-abundance-bubble.png", plot=accum_rel, width=20, height=20, units=c("cm"))
ggsave(filename="figs/ebpr-time-series-bubble.png", plot=fullset, width=40, height=20, units=c("cm"))
