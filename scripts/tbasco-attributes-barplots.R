library(tidyverse)
library(RColorBrewer)

# counts for all genomes 
attribute_counts <- read.csv("results/tbasco_results/all_genome_counts.csv")

attribute_counts$General_category <- gsub("#N/A", "Other Categories", attribute_counts$General_category)

attribute_info <- attribute_counts %>% 
  select(Trait_Attribute, General_category, Specific_category, Module_description)

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
  select(Trait_Attribute, Trait, Attribute, Genomes, filename, General_category, Specific_category, Module_description) %>% 
  filter(Genomes > 0)

lineage_table[is.na(lineage_table)] <- "Other Categories"

lineage.labs <- c("Actinobacteria (n=10)", "Alphaproteobacteria (n=15)", "Bacteroidetes (n=9)", "Betaproteobacteria (n=9)", "Gammaproteobacteria (n=4)")
names(lineage.labs) <- c("actinobacteria", "alphaproteobacteria", "bacteroidetes", "betaproteobacteria", "gammaproteobacteria")

lineage_attributes_plot <- lineage_table %>% ggplot(aes(x=factor(Genomes), fill=General_category)) + geom_bar() + facet_wrap(~ filename, scales=c("free_x"), labeller = labeller(filename=lineage.labs), nrow=1) + scale_y_continuous(expand=c(0,0), limits=c(0, 700), breaks=seq(0,700,100)) + scale_fill_brewer(palette="Set3") + labs(fill = "KEGG Module Category") + xlab("# of Genomes with Attribute") + ylab("Count") + theme_bw() + theme(legend.position="bottom")

ggsave("figs/all-genomes-attributes-barplot.png", all_genomes_attribute_barplot, width=13, height=6, units=c("in"))
ggsave("figs/lineage-attributes-barplot.png", lineage_attributes_plot, width=13, height=4, units=c("in"))

brewer.pal(10, "Set3")

#####################################
# Analysis of modules that are in core vs niche differentiating groups
#####################################

# Core >19 genomes 

attribute_counts %>% 
  filter(Genomes > 19) %>% 
  select(General_category) %>% 
  count(General_category, sort=TRUE)

attribute_counts %>% 
  filter(Genomes > 19) %>% 
  filter(General_category !=c("Other categories")) %>% 
  select(Specific_category) %>% 
  count(Specific_category, sort=TRUE)

core_attribute_modules <- attribute_counts %>% 
  filter(Genomes > 19) %>% 
  filter(General_category !=c("Other categories")) %>% 
  select(Module_description) %>% 
  count(Module_description, sort=TRUE)

attribute_counts %>% 
  filter(Genomes > 19) %>% 
  select(Trait) %>% 
  count(Trait, sort=TRUE)

# Niche differentiating, less than 10? 

attribute_counts %>% 
  filter(Genomes < 10) %>% 
  filter(General_category !=c("Other Categories")) %>% 
  select(General_category) %>% 
  count(General_category, sort=TRUE)

attribute_counts %>% 
  filter(Genomes < 10) %>% 
  filter(General_category !=c("Other Categories")) %>% 
  select(Specific_category) %>% 
  count(Specific_category, sort=TRUE)

differentiating_attribute_modules <- attribute_counts %>% 
  filter(Genomes < 10) %>% 
  filter(General_category !=c("Other Categories")) %>% 
  select(Module_description) %>% 
  count(Module_description, sort=TRUE) %>% 
  filter(n > 10)

# write out core and ND tables for all 

write_tsv(core_attribute_modules, "results/tbasco_results/core-attributes-summary.tsv")
write_tsv(differentiating_attribute_modules, "results/tbasco_results/differentiating-attributes-summary.tsv")

# core among different lineages 
colnames(lineage_table)[5] <- c("lineage")
# actinobacteria
core_actino <- lineage_table %>% 
  filter(lineage == c("actinobacteria")) %>% 
  filter(Genomes > 6) %>% 
  filter(General_category != c("Other Categories")) %>% 
  select(Module_description) %>% 
  count(Module_description, sort=TRUE) 

nd_actino <- lineage_table %>% 
  filter(lineage == c("actinobacteria")) %>% 
  filter(Genomes < 6) %>% 
  filter(General_category != c("Other Categories")) %>% 
  select(Module_description) %>% 
  count(Module_description, sort=TRUE) 

# alphaproteobacteria 

core_alpha <- lineage_table %>% 
  filter(lineage == c("alphaproteobacteria")) %>%
  filter(Genomes > 10) %>% 
  filter(General_category != c("Other categories")) %>% 
  select(Module_description) %>% 
  count(Module_description, sort=TRUE) 

nd_alpha <- lineage_table %>% 
  filter(lineage == c("alphaproteobacteria")) %>% 
  filter(Genomes < 10) %>% 
  filter(General_category != c("Other categories")) %>% 
  select(Module_description) %>% 
  count(Module_description, sort=TRUE) 

# bacteroidetes
core_bacteroidetes <- lineage_table %>% 
  filter(lineage == c("bacteroidetes")) %>%
  filter(Genomes > 6) %>% 
  filter(General_category != c("Other categories")) %>% 
  select(Module_description) %>% 
  count(Module_description, sort=TRUE) 

nd_bacteroidetes <- lineage_table %>% 
  filter(lineage == c("bacteroidetes")) %>% 
  filter(Genomes < 6) %>% 
  filter(General_category != c("Other categories")) %>% 
  select(Module_description) %>% 
  count(Module_description, sort=TRUE) 

# betaproteobacteria

core_beta <- lineage_table %>% 
  filter(lineage == c("betaproteobacteria")) %>% 
  filter(Genomes > 6) %>% 
  filter(General_category != c("Other categories")) %>% 
  select(Module_description) %>% 
  count(Module_description, sort=TRUE) 

nd_beta <- lineage_table %>% 
  filter(lineage == c("betaproteobacteria")) %>% 
  filter(Genomes < 6) %>% 
  filter(General_category != c("Other categories")) %>% 
  select(Module_description) %>% 
  count(Module_description, sort=TRUE) 

# write out lineage core and ND tables 
write_tsv(core_alpha, "results/tbasco_results/core-alpha-attributes.tsv")
write_tsv(nd_alpha, "results/tbasco_results/nd-alpha-attributes.tsv")
write_tsv(core_actino, "results/tbasco_results/core-actino-attributes.tsv")
write_tsv(nd_actino, "results/tbasco_results/nd-actino-attributes.tsv")
write_tsv(core_bacteroidetes, "results/tbasco_results/core-bact-attributes.tsv")
write_tsv(nd_bacteroidetes, "results/tbasco_results/nd-bact-atributes.tsv")
write_tsv(core_beta, "results/tbasco_results/core-beta-attributes.tsv")
write_tsv(nd_beta, "results/tbasco_results/nd-beta-attributes.tsv")
