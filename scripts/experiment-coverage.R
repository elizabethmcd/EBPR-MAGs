library(tidyverse)

coverage_data <- read_tsv("results/2013_binning/mapping/all_covg_results.txt", col_names = FALSE)
colnames(coverage_data) <- c("bin", "sample", "coverage")

metadata <- read.csv("results/R1R2-EBPR-MAGs-table.csv") %>% 
  select(Bin, Code)
colnames(metadata) <- c("bin", "code")

aggregate(coverage ~ bin, coverage_data, mean)

experiment_coverage <- coverage_data %>% 
  filter(sample == '2013-05-28-EBPR') %>% 
  select(bin, coverage)

coverage_metadata <- left_join(experiment_coverage, metadata) %>% 
  filter(code != "CAPIA" & code !="CAPIIA")

write.csv(coverage_metadata, "results/2013_binning/mapping/flanking-bins-experiment-coverage.csv", quote = FALSE, row.names = FALSE)
