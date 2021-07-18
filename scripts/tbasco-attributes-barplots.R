library(tidyverse)
library(RColorBrewer)

# counts for all genomes 
attribute_counts <- read.csv("results/tbasco_results/all_genome_counts.csv")

attribute_counts$General_category <- gsub("#N/A", "Other Categories", attribute_counts$General_category)

attribute_info <- attribute_counts %>% 
  select(Trait_Attribute, General_category)

all_genomes_attribute_barplot <- attribute_counts %>% ggplot(aes(x=factor(Genomes), fill=General_category)) + geom_bar() + scale_y_continuous(expand=c(0,0), limits=c(0,200), breaks=seq(0,200,25)) + scale_fill_brewer(palette="Set3") + labs(fill = "KEGG Module Category") + xlab("# of Genomes with Attribute") + ylab("Count") + theme_bw() + theme(legend.position="bottom")
all_genomes_attribute_barplot

# counts within lineages 

files_path <- "results/tbasco_results"
attribute_files <- dir(files_path, pattern="_attributes.csv")

lineage_attributes <- data_frame(filename = attribute_files) %>% 
  mutate(file_contents = map(filename, ~ read_csv(file.path(files_path, .)))) %>% 
  unnest()

lineage_attributes$filename <- gsub("_attributes.csv", "", lineage_attributes$filename)

lineage_table <- left_join(lineage_attributes, attribute_info) %>% 
  select(Trait_Attribute, Trait, Attribute, Genomes, filename, General_category) %>% 
  filter(Genomes > 0)

lineage_table[is.na(lineage_table)] <- "Other Categories"

lineage.labs <- c("Actinobacteria (n=10)", "Alphaproteobacteria (n=15)", "Bacteroidetes (n=9)", "Betaproteobacteria (n=9)", "Gammaproteobacteria (n=4)")
names(lineage.labs) <- c("actinobacteria", "alphaproteobacteria", "bacteroidetes", "betaproteobacteria", "gammaproteobacteria")

lineage_attributes_plot <- lineage_table %>% ggplot(aes(x=factor(Genomes), fill=General_category)) + geom_bar() + facet_wrap(~ filename, scales=c("free_x"), labeller = labeller(filename=lineage.labs), nrow=1) + scale_y_continuous(expand=c(0,0), limits=c(0, 700), breaks=seq(0,700,100)) + scale_fill_brewer(palette="Set3") + labs(fill = "KEGG Module Category") + xlab("# of Genomes with Attribute") + ylab("Count") + theme_bw() + theme(legend.position="bottom")

ggsave("figs/all-genomes-attributes-barplot.png", all_genomes_attribute_barplot, width=13, height=6, units=c("in"))
ggsave("figs/lineage-attributes-barplot.png", lineage_attributes_plot, width=13, height=4, units=c("in"))

brewer.pal(10, "Set3")
