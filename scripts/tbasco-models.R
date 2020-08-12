library(TbasCO)
library(tidyverse)

model <- '3300009517-bin.35'
Module_Names = names(ret$RNAseq.data$features$annotation.db$module.dict)
Model_Module_List_All <- Model_Module(RNAseq.data = ret$RNAseq.data, trait.attributes = ret$trait.attributes, model, bkgd.traits = bkgd.traits, Module_Names = Module_Names)

plot(sort(Model_Module_List_All$Fish_Backgrounds_trimmed))
margins =c(2,3)
Plot_Model_Module(Model_Module_List_All, model, Module_Names,margins, sortbygenome='3300009517-bin.35')

runella_comparisons <- as.data.frame(Model_Module_List_All[["Fish_Backgrounds_trimmed"]])

write.csv(runella_comparisons, "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R1R2/EBPR-MAGs/results/2013_binning/annotations/runella/runella_tbasco_module_comparisons.csv", quote=FALSE)
write.csv(test2, "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R1R2/EBPR-MAGs/results/2013_binning/annotations/runella/runella_tbasco_pairwise_comparisons.csv", quote=FALSE)

colnames(runella_comparisons) <- c("z_score")

runella_comparisons <- rownames_to_column(runella_comparisons, "module")

##########################################
# Accumulibacter pairwise comparisons 

ac_model <- '3300026302-bin.3'
Module_Names = names(ret$RNAseq.data$features$annotation.db$module.dict)
Ac_Model_Module_List_All <- Model_Module(RNAseq.data = ret$RNAseq.data, trait.attributes = ret$trait.attributes, ac_model, bkgd.traits = bkgd.traits, Module_Names = Module_Names)

ac_pairwise <- as.data.frame(Ac_Model_Module_List_All[["Model_Comparison_Matrix"]])
write.csv(ac_pairwise, "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R1R2/EBPR-MAGs/results/2013_binning/annotations/CAPIIA_tbasco_pairwise_comparisons.csv", quote=FALSE)

ac_model <- rownames_to_column(ac_pairwise, "bin")
ac_bin <- ac_model %>% 
  filter(bin == "3300026286-bin.31") %>% 
  select(-bin)

ac_long <- gather(ac_bin, key="module", value="z_score")

ac_clean <- ac_long %>% 
  drop_na() %>% 
  filter(z_score > 1 | z_score < -1)

ac_high_sig <- ac_clean %>% 
  filter(z_score > 2 | z_score < -2)

#################################
# Merge lists with KEGG module descriptions list
kegg_modules <- read.delim("/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R1R2/EBPR-MAGs/metadata/kegg_modules_list.tsv", sep="\t", header=FALSE)
colnames(kegg_modules) <- c("module_name", "description")

ac_high_clean <- ac_high_sig %>% 
  separate(module, into=c("module_name"), sep="\\.")
ac_high_sig_names <- left_join(ac_high_clean, kegg_modules)

ac_sig_clean <- ac_clean %>% 
  separate(module, into=c("module_name"), sep="\\.")
ac_all_sig_names <- left_join(ac_sig_clean, kegg_modules)

ac_med_sig <- ac_all_sig_names %>% 
  filter(z_score < -1.5 | z_score > 1.5)

write.csv(ac_all_sig_names, "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R1R2/EBPR-MAGs/results/2013_binning/annotations/accumulibacter/CAP_all_sig_comparisons.csv", row.names = FALSE, quote=FALSE)
write.csv(ac_med_sig, "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R1R2/EBPR-MAGs/results/2013_binning/annotations/accumulibacter/CAP_med_sig_comparisons.csv", row.names = FALSE, quote=FALSE)
write.csv(ac_high_sig_names, "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R1R2/EBPR-MAGs/results/2013_binning/annotations/accumulibacter/CAP_high_sig_comparisons.csv", row.names = FALSE, quote=FALSE)
