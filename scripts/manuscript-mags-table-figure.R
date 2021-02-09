# creating supplementary table of 16 high-coverage genomes

library(formattable)
library(dplyr)

# stats of full dataset
stats <- read.csv("/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/results/stats/refined-bin-stats.txt", sep="\t", header=TRUE)

# classifications of 16 high coverage MAGs 
classf <- read.csv("/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/results/stats/high-coverage-mags-gtdb-ncbi-classifications.csv", sep=",", header=TRUE)

# merge to get stats and names of only high covg bins
full <- left_join(classf, stats)

# size in Mbp
full$size_Mbp <- full$size / 1000000
full_no_qual <- full[,c(-6,-9,-11,-12)]

# sig figs
full_no_qual[6:7] <- round(full_no_qual[6:7], digits=2)
full_no_qual[9] <- round(full_no_qual[9], digits=2)

# reorder 
final <- full_no_qual[,c(2,6,7,9,3,4,8)]

# export
write.csv(final, file="~/Desktop/ebpr-mags-table.csv")
