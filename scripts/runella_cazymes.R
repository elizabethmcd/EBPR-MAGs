# runella CAZymes 
library(tidyverse)
library(reshape2)
library(viridis)

# counts of cazymes
runella_counts <- read.csv("results/2013_binning/annotations/runella/runella-cazymes-counts-no-annotations.csv", header=FALSE)
colnames(runella_counts) <- c("locus_tag", "an1", "an2", "an3", "aer1", "aer2", "aer3")

# cazyme annotations
runella_cazymes <- read.delim("results/2013_binning/annotations/runella/runella_cazymes_loci.txt", header=FALSE)
colnames(runella_cazymes) <- c("locus_tag", "cazyme_ID")

# merged
runella_cazyme_table <- left_join(runella_counts, runella_cazymes)

# heatmap of expression
runella_cazyme_table$total = rowSums(runella_cazyme_table[2:7])

top_loci <- runella_cazyme_table %>% 
  select(locus_tag, total) %>% 
  filter(total > 20) %>% 
  pull(locus_tag)

top_sets <- runella_cazyme_table %>% 
  filter(locus_tag %in% top_loci) %>% 
  select(locus_tag, an1, an2, an3, aer1, aer2, aer3)

cazy_melt <- melt(top_sets, id.vars="locus_tag", variable="condition")
ggplot(cazy_melt, aes(x=condition, y=locus_tag, fill=value)) + geom_tile(color="white") + scale_fill_viridis(option="viridis", alpha=1, begin=0, end=1, direction=-1)


top_heatmap <- ggplot(top_ts, aes(x=code, y=fct_rev(date), fill=value)) + geom_tile(color="white") + scale_fill_viridis(option="viridis", alpha=1, begin=0, end=1, direction=-1)
