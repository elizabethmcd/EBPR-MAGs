library(tximport)
library(tximportData)
library(readr)
library(tibble)
library(dplyr)
library(ggplot2)

# kallisto files
dir <- "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/results/transcriptomics/kallisto_results"
samples <- read.table(file.path(dir, "samples.txt"), header=TRUE)
files <- file.path(dir, samples$experiment, "abundance.h5")
names(files) <- paste0("sample", 1:6)
txi.kallisto <- tximport(files, type="kallisto", txOut = TRUE)
counts <- as.data.frame(txi.kallisto)
finalcounts <- rownames_to_column(counts, var="ID")
counttable <- finalcounts[, c(1:7)]
colnames(counttable) <- c("Locus_Tag", "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6")
write_delim(counttable, "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/results/transcriptomics/trans-mapping/ebpr-kallisto-tpm-counts-norm.tsv", delim="\t")

# merge counts and annotations

countsFile = read.csv("/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/results/transcriptomics/trans-mapping/2020-02-21-normalized-tmp-names.tsv", header=FALSE, sep="\t")
annotFile = read.csv("/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/results/annotations/2019-07-29-KO-redone/ebpr-kofamkoala-annots-sig-mod-nodups-ko-list.txt", header=FALSE, sep="\t")
colnames(countsFile) <- c("Bin", "Locus_Tag", "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6")
colnames(annotFile) <- c("Locus_Tag", "Annotation")
tpmTable = left_join(annotFile, countsFile)
write_delim(tpmTable, "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/results/transcriptomics/trans-mapping/2020-02-21-tpm-normalized-counts-annots.tsv", delim="\t")
countsTable <- rawTable[,c(2:8,9,1)]

# check to see merged correctly 
annotCounts <- countsTable %>% select(Bin, Annotation) %>% group_by(Bin) %>% mutate(Annotation.count = n()) %>% slice(1) %>% unique()

# high mapping bins

highList = c('3300009517-bin.1',
             '3300009517-bin.12',
             '3300009517-bin.13',
             '3300009517-bin.3',
             '3300009517-bin.31',
             '3300009517-bin.42',
             '3300009517-bin.47',
             '3300009517-bin.6',
             '3300026282-bin.4',
             '3300026284-bin.9',
             '3300026288-bin.43',
             '3300026302-bin.10',
             '3300026302-bin.32',
             '3300026302-bin.46',
             '3300026302-bin.62',
             '3300026303-bin.42')

# check sums and merged annotations
highCovg = countsTable %>% filter(Bin %in% highList)
HighannotCounts <- highCovg %>% select(Bin, Annotation) %>% group_by(Bin) %>% mutate(Annotation.count = n()) %>% slice(1) %>% unique()
countSums <- aggregate(highCovg[2:7], list(highCovg$Bin), sum)
countTotals <- cbind.data.frame(countSums$Group.1, rowSums(countSums[2:7]))
colnames(countTotals) <- c("Bin", "Total_Raw_Counts")
countTotals[order(countTotals$Total_Raw_Counts),]

# save files
write.table(countsTable, "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/raw-data/tbasco-tables/2019-08-10-full-tbasco-input-table-updated-annotations.tsv", sep=";", row.names=FALSE, quote=FALSE)
write.table(highCovg, "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/raw-data/tbasco-tables/2019-08-10-high-covg-tbasco-input-table-updated-annotations.tsv", sep=";", row.names=FALSE, quote=FALSE)

# RPKM normalization
# merge with Prokka annotations to normalize
prokka <- read.delim("/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/raw-data/annotations/ebpr-prokka-annotations.txt", sep="\t", header=FALSE)
colnames(prokka) <- c("Locus_Tag", "prokka_annotation", "size_bp", "accession")
counts_all_annots <- left_join(countsTable, prokka)
counts_annots_high_covg <- left_join(highCovg, prokka)

# remove rows with gene sizes as NA to throw out rRNA/tRNA hits that aren't annotated and throw off analyses
counts_no_nas <- counts_annots_high_covg[!is.na(counts_annots_high_covg$size_bp), ]

# divide by gene lengths
counts_no_nas$size_kbp <- counts_no_nas$size_bp / 1000
counts_rpk = as.data.frame(lapply(counts_no_nas[,c(2:7)], function(x) {
  (x / counts_no_nas$size_kbp)
  }))

# per million factor
counts_rpk["PM", ] = (colSums(counts_rpk[1:6]) / 100000)

# divide by per million factor
counts_tpm = as.data.frame(lapply(counts_rpk, function(x) x/tail(x, 1) ))
counts_tpm = counts_tpm[-63507,]

# put in genome names and locus tags
counts_tpm$genome_name <- counts_no_nas$Bin
counts_tpm$locus_tag <- counts_no_nas$Locus_Tag
counts_totals <- aggregate(counts_tpm[1:6], list(counts_tpm$genome_name), sum)

# bin metadata
metadata <- read.csv("/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/results/mapping/ebpr-rank-abundance-classifications.csv")
bin_names <- metadata %>% select("bin", "Highest_Classf")
colnames(counts_totals)[1] = "bin"
counts_names <- left_join(counts_totals, bin_names)

# expression averages
expression_averages <- cbind.data.frame(counts_names$bin, (rowSums(counts_names[2:7]) / 6), counts_names$Highest_Classf)
colnames(expression_averages) <- c("Genome", "Average_Expression", "Classification")

# plot
expression_averages %>% ggplot(aes(x=reorder(Genome,-Average_Expression), y=Average_Expression, fill=Classification)) + geom_col() + labs(x="Genome", y="Total Expression") + theme_classic()

# merge TPM normalization with KO annotations
colnames(counts_tpm)[8] <- "Locus_Tag"
tpm_kos <- left_join(counts_tpm, annotFile)
tpm_kos[1:6] <- round(tpm_kos[1:6], digits=3)

# fix column names
colnames(tpm_kos) <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6", "Bin", "Locus_Tag", "Annotation")
tpm_kos <- tpm_kos[, c(8,1:6,9,7)]

# export TPM normalized dataframe w/ KO annotations
write.table(tpm_kos, "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/raw-data/tbasco-tables/2019-08-10-high-covg-tpm-normalized-fixed-annotations.tsv", sep=";", row.names=FALSE, quote=FALSE)


##################################
# Get TPM normalization from kallisto files directly
# Merge with bin names
# Average total expression across 6 samples per genome

tax = read.csv("/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/results/stats/ebpr-mags-phyla-classf-stats.txt", sep="\t", header=FALSE)
colnames(tax) = c("Genome", "Group")
tpm = finalcounts[1:7]
write_delim(tpm, "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/results/transcriptomics/trans-mapping/ebpr-kallisto-tpm-normalized-counts.tsv", delim="\t")
tpm_names = read.csv("/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/results/transcriptomics/trans-mapping/2020-02-19-tpm-norm-names.tsv", sep="\t" )
colnames(tpm_names) = c("Genome", "Locus_Tag", "sample1", "sample2", "sample3", "sample4", "sample5", "sample6")
counts = tpm_names[-2]
countSums <- aggregate(counts[2:7], list(counts$Genome), sum)
countAvg <- cbind.data.frame(countSums$Group.1, (rowSums(countSums[2:7]) / 6))
colnames(countAvg) = c("Genome", "AvgExp")
countNames = left_join(tax, countAvg)
# relative expression plot
relexpplot = ggplot(countNames, aes(x=reorder(Genome,-AvgExp), y=AvgExp, fill=Group)) + geom_col() + labs(x="Genome", y="Avrage Expression Across Time-Series") + theme_classic() + scale_fill_brewer(palette="Paired") 
relexpplot
ggsave(filename="/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/figs/avg-rel-exp-all-bins.png", plot=relexpplot, width=40, height=20, units=c("cm"))

# looking at total raw counts across all 6 samples
raw = read.csv("/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/results/transcriptomics/trans-mapping/ebpr-raw-counts-names.tsv", header=FALSE, sep="\t")
colnames(raw) = c("Genome", "Locus_Tag", "sample1", "sample2", "sample3", "sample4", "sample5", "sample6")
raw_counts = raw[-2]
rawSums = aggregate(raw_counts[2:7], list(raw_counts$Genome), sum)
allSums = cbind.data.frame(rawSums$Group.1, (rowSums(rawSums[2:7])))
colnames(allSums) = c("Genome", "Total_Counts")
allSums_names= left_join(tax, allSums)
write_csv(allSums_names, "~/Desktop/counts_results.csv", quote=FALSE)
