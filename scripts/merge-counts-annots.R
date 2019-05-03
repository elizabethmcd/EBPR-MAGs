library(dplyr)

# merge counts and annotations

countsFile = read.csv("raw-data/ebpr-raw-counts-names.tsv", header=FALSE, sep="\t")
annotFile = read.csv("raw-data/ebpr-kofamAnnotations.txt", header=FALSE, sep="\t")
colnames(countsFile) <- c("Bin", "Locus_Tag", "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6")
colnames(annotFile) <- c("Locus_Tag", "Annotation")
rawTable <- left_join(countsFile, annotFile)
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
highCovg = countsTable %>% filter(Bin %in% highList)
HighannotCounts <- highCovg %>% select(Bin, Annotation) %>% group_by(Bin) %>% mutate(Annotation.count = n()) %>% slice(1) %>% unique()

# save files
write.table(countsTable, "raw-data/full-tbasco-input-table.tsv", sep=";", row.names=FALSE, quote=FALSE)
write.table(highCovg, "raw-data/high-covg-tbasco-input-table.tsv", sep=";", row.names=FALSE, quote=FALSE)