library(tidyverse)

counts = read.csv('results/2013_transcriptomes/2013_R1R2_ebpr_raw_counts_kallisto.csv')
norm = read.csv('results/2013_transcriptomes/2013_R1R2_ebpr_normalized_TPM_counts_kallisto.csv')
annotations = read.delim('results/2013_binning/annotations/R1R2-formatted-annotations.txt', header=FALSE)
colnames(annotations) = c('Bin', 'Locus_Tag', 'KO')
merged = left_join(counts,annotations) %>% select(Locus_Tag, B_15min_Anaerobic, D_52min_Anaerobic, F_92min_Anaerobic, H_11min_Aerobic, J_51min_Aerobic, N_134min_Aerobic, KO, Bin)
colnames(merged)[8] = c('Annotation')
norm_merged = left_join(norm,annotations) %>% select(Locus_Tag, B_15min_Anaerobic, D_52min_Anaerobic, F_92min_Anaerobic, H_11min_Aerobic, J_51min_Aerobic, N_134min_Aerobic, KO, Bin)
colnames(norm_merged)[8] = c('Annotation')

write_delim(merged, 'results/2013_transcriptomes/tables/2013_R1R2_raw_table.csv', delim = ';')
write_delim(norm_merged, 'results/2013_transcriptomes/tables/2013_R1R2_norm_table.csv', delim=';')
