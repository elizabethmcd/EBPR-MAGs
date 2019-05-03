library(dplyr)

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  counts <- args[1]
  annotation <- args[2]
  output <- args[3]
}


mergeDatasets <- function(counts, annotation) {
  countsFile <- read.csv(counts, header=TRUE, sep="\t")
  annotFile <- read.csv(annotation, header=FALSE, sep="\t")
  colnames(annotFile) <- c("Bin", "Locus_Tag", "Annotation")
  rawTable <- left_join(countsFile, annotFile)
  countsTable <- rawTable[,c(1:7,9,8)]
  countsOfAnnotations <- aggregate(countsTable[8], list(countsTable$Annotation))
}

countsFile = read.csv("raw-data/ebpr-raw-counts-names.tsv", header=FALSE, sep="\t")
annotFile = read.csv("raw-data/ebpr-kofamAnnotations.txt", header=FALSE, sep="\t")
colnames(countsFile) <- c("Bin", "Locus_Tag", "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6")
colnames(annotFile) <- c("Locus_Tag", "Annotation")
rawTable <- left_join(countsFile, annotFile)
countsTable <- rawTable[,c(2:8,9,1)]

annotCounts <- rawTable %>% select(Bin, Annotation) %>% group_by(Bin) %>% mutate(Annotation.count = n()) %>% slice(1) %>% unique()
naLines <- rawTable %>% filter(Bin == "")
countsTable %>% select("Annotation") %>% unique %>% count()
write.table(countsTable, "raw-data/2019-04-26-input-table.tsv", sep=";", row.names=FALSE, quote=FALSE)




main()
mergeDatasets()