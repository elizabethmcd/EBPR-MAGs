library(tidyverse)

mapping = read_delim('results/2013_mapping/R1R2-abundance.txt', delim="\t", col_names=FALSE)
stats = read_delim('results/2013_binning/R1R2-May2013-final-rep-species-bins-stats-table.csv', delim=',')
colnames(mapping) = c('Bin', '2013-05-13', '2013-05-23', '2013-05-28')
merged = left_join(stats, mapping)
write.csv(merged, 'results/2013_binning/2013_R1R2_finalBins_stats_mapping.csv', quote=FALSE, row.names = FALSE)

# manual classfs added
ebpr = read.csv('results/2013_R1R2_finalBins_stats_mapping_classfs.csv')

# TPM transcriptional results added
counts = read.csv('results/2013_transcriptomes/2013_R1R2_ebpr_normalized_TPM_counts_kallisto.csv')

# quality figure
qual = ebpr %>% select(Highest_Classf, comp, contam)
qualFig = qual %>% ggplot(aes(x=comp, y=contam, color=Highest_Classf)) + geom_point(size=4) + scale_color_brewer(palette="Paired") + theme_classic()
qualFig

# averaging relative abundance across all may samples and 28th may day for when transcriptional experiment was done
abundance = ebpr %>% select(Bin, Highest_Classf, X5.13.13, X5.23.13, X5.28.13)
abundance$average = (abundance$X5.13.13 + abundance$X5.23.13 + abundance$X5.28.13) / 3
may2013 = abundance %>% select(Bin, Highest_Classf, X5.28.13)

# figures
avgFig = abundance %>% ggplot(aes(x=reorder(Bin, -average), y=average, fill=Highest_Classf)) + geom_col() + labs(x="Genome", y="% Average Relative Abundance") + theme_classic() + scale_fill_brewer(palette="Paired")
avgFig
mayFig = may2013 %>% ggplot(aes(x=reorder(Bin, -X5.28.13), y=X5.28.13, fill=Highest_Classf)) + geom_col() + labs(x="Genome", y="% Relative Abundance") + theme_classic() + scale_fill_brewer(palette="Paired")
mayFig

# parsing transcriptional results
sumCounts = aggregate(counts[3:8], list(counts$Bin), sum)
sumCounts$Avg = rowMeans(sumCounts[,-1])
colnames(sumCounts)[1] = "Bin"

# merge with classf names
names = ebpr %>% select(Bin, Highest_Classf)
countsClassf = left_join(names, sumCounts)
avgCounts = countsClassf %>% select(Bin, Highest_Classf, Avg)
countsFig = avgCounts %>% ggplot(aes(x=reorder(Bin, -Avg), y=Avg, fill=Highest_Classf)) + geom_col() + labs(x="Genome", y="Average Relative Activity in Transcripts per Million [TPM]") + theme_classic() + scale_fill_brewer(palette="Paired")


  