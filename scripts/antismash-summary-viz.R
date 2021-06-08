library(tidyverse)

# EBPR MAGs antismash summaries

secmet <- read.csv("results/antismash_annotations/EBPR-MAGS-antismash-summs-results.csv")
secmet[is.na(secmet)] <- 0

counts <- secmet %>% select(-Bin, -Code, -Classification, -Total) %>% 
  pivot_longer(-Phylum, names_to="type", values_to="count")

counts_filtered <- counts %>% filter(Phylum != "Patescibacteria") %>% filter(Phylum != "Firmicutes")

secmet_plot <- counts_filtered %>% ggplot(aes(x=fct_rev(Phylum), y=count, fill=type)) + geom_bar(stat="identity") + scale_fill_brewer(palette="Spectral", labels=c("Acyl Amino Acid (5)", "Bacteriocin (32)", "Betalactone (14)", "Lanthipeptide (3)", "Lassopeptide (8)", "NRPS/PKS (60)", "Other (LAP, Aryl Polyene, Homserine Lactone) (33)", "Terpene (48)"), name="Secondary Metabolite Type") + coord_flip() + scale_y_continuous(expand=c(0,0), limits=c(0,115), breaks=seq(0,115,10)) + scale_x_discrete(labels=c("Proteobacteria (34)", "Melainabacteria (1)", "Gemmatimonadetes (1)", "Chloroflexi (2)", "Bacteroidetes (13)", "Armatimonadetes (1)", "Actinobacteria (12)")) + theme_classic() +  ggtitle("Secondary Metabolites Predicted Among EBPR SBR Genomes") + ylab("\nNumber of Clusters\n") + xlab("\nGenome Phylum\n") + theme(axis.text.y=element_text(size=12, face="italic"), axis.text.x=element_text(size=12), legend.position=c(0.6, 0.5), plot.title=element_text(size=16, face="bold"), legend.title=element_text(size=14, face="bold"), legend.text=element_text(size=12), axis.title.x=element_text(size=14, face="bold"), axis.title.y=element_text(size=14, face="bold"))

ggsave(filename="figs/EBPR-SecMets-Summs.png", plot=secmet_plot, width=15, height=6, units=c("in"))
