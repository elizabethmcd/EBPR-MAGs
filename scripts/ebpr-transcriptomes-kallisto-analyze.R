library(tximport)
library(readr)
library(tibble)
library(dplyr)
library(ggplot2)

# kallisto files
dir <- "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R1R2/EBPR-MAGs/results/2013_transcriptomes/"
samples <- read.table(file.path(dir, "samples.txt"), header=TRUE)
files <- file.path(dir, samples$experiment, "abundance.h5")
names(files) <- paste0("sample", 1:6)
txi.kallisto <- tximport(files, type="kallisto", txOut = TRUE)
counts <- as.data.frame(txi.kallisto)
finalcounts <- rownames_to_column(counts, var="ID")
write.csv(finalcounts, 'results/2013_transcriptomes/2013R1R2-ebpr-bins-kallisto-count-table.csv', row.names = FALSE, quote=FALSE)

# tables of raw and normalized counts from txi dataframe
# raw
rawcounts <- as.data.frame(txi.kallisto$counts)
rawTable <- rownames_to_column(rawcounts, var="ID")
rawTable.split <- rawTable %>% separate(ID, c("genome"), sep='_') %>% cbind(rawTable$ID)
colnames(rawTable.split) <- c('Bin', 'B_15min_Anaerobic', 'D_52min_Anaerobic', 'F_92min_Anaerobic', 'H_11min_Aerobic', 'J_51min_Aerobic', 'N_134min_Aerobic', 'Locus_Tag')
rawOut <- rawTable.split %>% select(Bin, Locus_Tag, B_15min_Anaerobic, D_52min_Anaerobic, F_92min_Anaerobic, H_11min_Aerobic, J_51min_Aerobic, N_134min_Aerobic)
write.csv(rawOut, 'results/2013_transcriptomes/2013_R1R2_ebpr_raw_counts_kallisto.csv', quote = FALSE, row.names = FALSE)

# normalized
normcounts <- as.data.frame(txi.kallisto$abundance)
normTable <- rownames_to_column(normcounts, var='ID')
norm.split <- normTable %>% separate(ID, c("genome"), sep='_') %>% cbind(normTable$ID)
colnames(norm.split) <- c('Bin', 'B_15min_Anaerobic', 'D_52min_Anaerobic', 'F_92min_Anaerobic', 'H_11min_Aerobic', 'J_51min_Aerobic', 'N_134min_Aerobic', 'Locus_Tag')
normOut <- norm.split %>% select(Bin, Locus_Tag, B_15min_Anaerobic, D_52min_Anaerobic, F_92min_Anaerobic, H_11min_Aerobic, J_51min_Aerobic, N_134min_Aerobic)
write.csv(normOut, 'results/2013_transcriptomes/2013_R1R2_ebpr_normalized_TPM_counts_kallisto.csv', quote=FALSE, row.names = FALSE)

