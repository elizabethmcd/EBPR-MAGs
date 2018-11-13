# Analyzing metatranscriptomic experiments and EBPR MAGs expression 

library(tidyverse)

# Read in read count files for each experiment 
Ban = read.delim("raw-data/B_15min_Anaerobic-raw-counts.txt", sep="\t", header=FALSE, col.names = c("bin_contig", "locus_tag", "raw_counts"))
Dan = read.delim("raw-data/D_52min_Anaerobic-raw-counts.txt", sep ="\t", header=FALSE, col.names = c("bin_contig", "locus_tag", "raw_counts"))
Fan = read.delim("raw-data/F_92min_Anaerobic-raw-counts.txt", sep="\t", header=FALSE, col.names = c("bin_contig", "locus_tag", "raw_counts"))
Hae = read.delim("raw-data/H_11min_Aerobic-raw-counts.txt", sep="\t", header=FALSE, col.names = c("bin_contig", "locus_tag", "raw_counts"))
Jae = read.delim("raw-data/J_51min_Aerobic-raw-counts.txt", sep="\t", header=FALSE, col.names = c("bin_contig", "locus_tag", "raw_counts"))
Nae = read.delim("raw-data/N_134min_Aerobic-raw-counts.txt", sep="\t", header=FALSE, col.names = c("bin_contig", "locus_tag", "raw_counts"))

# Taxonomy file
taxonomy = read.csv("results/ebpr-rank-abundance-classifications.csv", header=TRUE)
taxNames = subset(taxonomy, select=c(Bin, Highest_Classf))

# Genome stats
genomeStats = read.delim("results/refined-bins-derep-checkm-stats.txt", sep=" ", header=FALSE, col.names=c("Bin", "Classification", "Size", "Completeness", "Redundancy"))
genomeStats = subset(genomeStats, select=c(Bin, Size))

Ban = separate(Ban, bin_contig, into=c("Bin", "Contig"), sep="_")
Dan = separate(Dan, bin_contig, into=c("Bin", "Contig"), sep="_")
Fan = separate(Fan, bin_contig, into=c("Bin", "Contig"), sep="_")
Hae = separate(Hae, bin_contig, into=c("Bin", "Contig"), sep="_")
Jae = separate(Jae, bin_contig, into=c("Bin", "Contig"), sep="_")
Nae = separate(Nae, bin_contig, into=c("Bin", "Contig"), sep="_")

# Example with 1st anaerobic stage sample 
Ban = subset(Ban, select = -c(Contig))
BanNames = left_join(Ban, taxNames)
BanNames$normalized_total = BanNames$raw_counts / 256007994        
BanGenomes = left_join(BanNames, genomeStats)
BanFull = transform(BanGenomes, normalized_total_genome = normalized_total / Size)
BanCounts = as.data.frame(spread(BanFull, key="Bin", value="normalized_total_genome", fill=0))
BanCounts = subset(BanCounts, select = -c(locus_tag, Highest_Classf, normalized_total, Size, raw_counts))
totalBanCounts = as.data.frame(lapply(BanCounts, as.numeric))
binNames = colnames(BanCounts)
colnames(totalBanCounts) = binNames
fullCounts = as.data.frame(colSums(totalBanCounts))
fullCounts = rownames_to_column(fullCounts, "Bin")
colnames(fullCounts) = c("Bin", "normalized_expression")
CountNames = left_join(fullCounts, taxNames)

# plot rank expression B
BPLOT <- ggplot(CountNames, aes(x=reorder(Bin,-normalized_expression), y=normalized_expression, fill=Highest_Classf)) + geom_col() + labs(x="Genome", y="Total Expression") + theme_classic()


# D timepoint
Dan = subset(Dan, select = -c(Contig))
DanNames = left_join(Dan, taxNames)
DanNames$normalized_total = DanNames$raw_counts / 202317784       
DanGenomes = left_join(DanNames, genomeStats)
DanFull = transform(DanGenomes, normalized_total_genome = normalized_total / Size)
DanCounts = as.data.frame(spread(DanFull, key="Bin", value="normalized_total_genome", fill=0))
DanCounts = subset(DanCounts, select = -c(locus_tag, Highest_Classf, normalized_total, Size, raw_counts))
totalDanCounts = as.data.frame(lapply(DanCounts, as.numeric))
DbinNames = colnames(DanCounts)
colnames(totalDanCounts) = DbinNames
fullCountsD = rownames_to_column(totalDanCounts, "Bin")
colnames(fullCountsD) = c("Bin", "normalized_expression")
CountNamesD = left_join(fullCountsD, taxNames)

# plot rank expression D
ggplot(CountNamesD, aes(x=reorder(Bin,-normalized_expression), y=normalized_expression, fill=Highest_Classf)) + geom_col() + labs(x="Genome", y="Total Expression") + theme_classic()

# F sample
Fan = subset(Fan, select = -c(Contig))
FanNames = left_join(Fan, taxNames)
FanNames$normalized_total = FanNames$raw_counts / 252399352       
FanGenomes = left_join(FanNames, genomeStats)
FanFull = transform(FanGenomes, normalized_total_genome = normalized_total / Size)
FanCounts = as.data.frame(spread(FanFull, key="Bin", value="normalized_total_genome", fill=0))
FanCounts = subset(FanCounts, select = -c(locus_tag, Highest_Classf, normalized_total, Size, raw_counts))
totalFanCounts = as.data.frame(lapply(FanCounts, as.numeric))
FbinNames = colnames(FanCounts)
colnames(totalFanCounts) = FbinNames
fullCountsF = rownames_to_column(totalFanCounts, "Bin")
colnames(fullCountsF) = c("Bin", "normalized_expression")
CountNamesF = left_join(fullCountsF, taxNames)

# plot rank expression F
ggplot(CountNamesF, aes(x=reorder(Bin,-normalized_expression), y=normalized_expression, fill=Highest_Classf)) + geom_col() + labs(x="Genome", y="Total Expression") + theme_classic()

# H sample 
Hae = subset(Hae, select = -c(Contig))
HanNames = left_join(Hae, taxNames)
HanNames$normalized_total = HanNames$raw_counts / 214039116       
HanGenomes = left_join(HanNames, genomeStats)
HanFull = transform(HanGenomes, normalized_total_genome = normalized_total / Size)
HanCounts = as.data.frame(spread(HanFull, key="Bin", value="normalized_total_genome", fill=0))
HanCounts = subset(HanCounts, select = -c(locus_tag, Highest_Classf, normalized_total, Size, raw_counts))
totalHanCounts = as.data.frame(lapply(HanCounts, as.numeric))
binNamesH = colnames(HanCounts)
colnames(totalHanCounts) = binNamesH
fullCountsH = as.data.frame(colSums(totalHanCounts))
fullCountsH = rownames_to_column(fullCountsH, "Bin")
colnames(fullCountsH) = c("Bin", "normalized_expression")
CountNamesH = left_join(fullCountsH, taxNames)

# plot rank expression H
HPLOT <- ggplot(CountNamesH, aes(x=reorder(Bin,-normalized_expression), y=normalized_expression, fill=Highest_Classf)) + geom_col() + labs(x="Genome", y="Total Expression") + theme_classic()

# N sample 
Nae = subset(Nae, select = -c(Contig))
NanNames = left_join(Nae, taxNames)
NanNames$normalized_total = NanNames$raw_counts / 279356366       
NanGenomes = left_join(NanNames, genomeStats)
NanFull = transform(NanGenomes, normalized_total_genome = normalized_total / Size)
NanCounts = as.data.frame(spread(NanFull, key="Bin", value="normalized_total_genome", fill=0))
NanCounts = subset(NanCounts, select = -c(locus_tag, Highest_Classf, normalized_total, Size, raw_counts))
totalNanCounts = as.data.frame(lapply(NanCounts, as.numeric))
binNamesN = colnames(NanCounts)
colnames(totalNanCounts) = binNamesN
fullCountsN = as.data.frame(colSums(totalNanCounts))
fullCountsN = rownames_to_column(fullCountsN, "Bin")
colnames(fullCountsN) = c("Bin", "normalized_expression")
CountNamesN = left_join(fullCountsN, taxNames)

# plot rank expression N
ggplot(CountNamesN, aes(x=reorder(Bin,-normalized_expression), y=normalized_expression, fill=Highest_Classf)) + geom_col() + labs(x="Genome", y="Total Expression") + theme_classic()

# save B and H plots, beginning of anaerobic and beginning of aerobic stages to investigate for now
ggsave(file="figs/B-EXPRESSION.png", BPLOT, width=40, height=20, units=c("cm"))
ggsave(file="figs/H-EXPRESSION.png", HPLOT, width=40, height=20, units=c("cm"))
