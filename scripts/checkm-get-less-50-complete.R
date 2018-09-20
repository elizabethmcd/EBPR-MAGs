# Visualiation of completion/contamination of genome bins

library(dplyr)
library(ggplot2)
library(data.table)
library(gtable)

# Read in concatenated lineage text file
lineages = read.table("results/EBPR-bins-checkm-results.txt", sep="", header=FALSE)
names(lineages) <- c("Bin", "Classification", "Genome_Length", "Completeness", "Redundancy")

# Set significant figures for completion and contamination 
lindf <- as.data.frame(lineages)
lindf$Genome_Length <- lindf$Genome_Length/1000000
head(lindf)

# Plot
p1 <- lindf %>% ggplot(aes(x=Completeness, y=Redundancy, color=Bin)) + geom_point(size=2)
p2 <- p1 + theme_bw() + theme(legend.position="none") 
# messy thing because lots of bins and decimal points that i'm not going to mess with right now

# Get a list of bins that are less than 50% complete to toss out of everything
less50 <- lindf %>% filter(Completeness < 50.00)
less50$Bin <- paste(less50$Bin, ".fna", sep="")
final <- less50$Bin
head(final)
write.table(final, file="less-than-50-complete-bins.txt", sep="\t", row.names=FALSE, quote=FALSE)

