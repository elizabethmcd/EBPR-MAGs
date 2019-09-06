library(TbasCO)

# File of raw counts and annotations
ebpr_counts = file.path("/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/EBPR-MAGs/results/transcriptomics/tbasco-tables/2019-08-10-high-covg-tbasco-input-table-updated-annotations.tsv")

# Normalization
normalization.features <- c(234232896, 183166236, 228746720, 198024002, 231567992, 259156166)

# RNA seq object
RNAseq.data <- Pre_process_input(ebpr_counts, annotation.db.path, normalize.method="simple", filter.method='MAD', filter.low.coverage = T, normalization.features)

# Classification of bins above cutoff
genome.taxonomy <- list('3300026303-bin.42'='f_Chitinophagaceae',
                        '3300026302-bin.62'='f_Obscuribacteraceae',
                        '3300026302-bin.47'='g_Pseudoxanthomonas',
                        '3300026302-bin.46'='g_Pseudoxanthomonas',
                        '3300026302-bin.32'='p_Bacteroidetes',
                        '3300026302-bin.20'='f_Burkholderiaceae',
                        '3300026302-bin.10'='f_Burkholderiaceae',
                        '3300026299-bin.22'='f_Dermatophilaceae',
                        '3300026288-bin.43'='g_Tetrasphaera',
                        '3300026288-bin.32'='f_Saprospiraceae',
                        '3300026284-bin.9'='Accumulibacter_IIA',
                        '3300026283-bin.28'='g_Sphingopyxis',
                        '3300026282-bin.4'='Accumulibacter_IA',
                        '3300009517-bin.7'='g_Salinibacterium',
                        '3300009517-bin.6'='g_Aquidulcibacter',
                        '3300009517-bin.47'='f_Rhodobacteraceae',
                        '3300009517-bin.42'='g_Rubrivivax',
                        '3300009517-bin.31'='g_Flavobacterium',
                        '3300009517-bin.30'='g_Chryseobacterium',
                        '3300009517-bin.3'='g_Tetrasphaera',
                        '3300009517-bin.13'='f_Burkholderiaceae',
                        '3300009517-bin.12'='g_Runella',
                        '3300009517-bin.1'='g_Leadbetterella')

# Distance metrics
PC <- function(rowA, rowB, RNAseq.features){
  return(cor(as.numeric(rowA[RNAseq.features$sample.columns]),
             as.numeric(rowB[RNAseq.features$sample.columns])
  )
  )
}
# Calculates the Normalized Rank Euclidean Distance
NRED <- function(rowA, rowB, RNAseq.features) {
  r.A <- as.numeric(rowA[ RNAseq.features$rank.columns ])
  r.B <- as.numeric(rowB[ RNAseq.features$rank.columns ])
  return(
    sum((r.A - r.B) * (r.A - r.B))
  )
}
# Combine multiple distance metrics to complement each other.
distance.metrics <- list("NRED" = NRED,
                         "PC"   = PC)
# the size of all modules in the library, usefull for calculating the background distribution
module_size_range <- lapply(RNAseq.data$features$annotation.db$module.dict,length) %>% as.numeric() %>% table() %>% names() %>% as.numeric()

# again, the size of all modules but using the maximum length of the disjunctive form. Slightly shorter
d_module_size_range_all<-c()
for (i in 1:length(sub_modules)) {
  d_module_size_range_all[i]<- max(sapply(sub_modules[[i]],length))
}
d_module_size_range<-sort(unique(d_module_size_range_all))

# Analysis
bkgd.individual <- Individual_Annotation_Background(RNAseq.data,N= 10000, metrics = distance.metrics, threads = 2)

bkgd.individual.Zscores <- Calc_Z_scores(bkgd.individual, distance.metrics)

bkgd.traits <- Random_Trait_Background(RNAseq.data, bkgd.individual.Zscores,N=10000,metrics=distance.metrics,threads=2)

pairwise.distances  <- Calc_Pairwise_Annotation_Distance(RNAseq.data, RNAseq.data$features$annotation.db, distance.metrics,bkgd.individual.Zscores, show.progress = T,threads = 4)

trait.attributes <- Identify_Trait_Attributes(RNAseq.data = RNAseq.data, pairwise.distances = pairwise.distances, annotation.db = RNAseq.data$features$annotation.db,threads = 2)

trait.attributes.pruned <- Prune_Trait_Attributes(trait.attributes, bkgd.traits, RNAseq.data, p.threshold = 0.05, pairwise.distances = pairwise.distances, bkgd.individual.Zscores = bkgd.individual.Zscores, annotation.db = RNAseq.data$features$annotation.db, trait_presence_absence = RNAseq.data$features$trait_presence_absence)

sbs.trait.attributes <- Traitattributes_To_Sbsmatrix(trait.attributes.pruned,RNAseq.data$features$bins)

# Plots of background distributions and metrics
Plot_Background_Modules(bkgd.traits)
Plot_Background_Individual_Genes(bkgd.individual.Zscores)
Plot_Metric_Comparison(bkgd.individual)
Plot_Redundancy_Traits(RNAseq.data)

# Redundancy figure
TnA_redundancy <- Calc_TnA_redundancy(RNAseq.data)
Most_redundant_order <- TnA_redundancy[order(TnA_redundancy[,3]),1]
stc <- string.to.colors(genome.taxonomy)
stc_highlights <- stc[match(rev(Most_redundant_order),names(stc))]
par(mfrow=c(4,5),mar=c(2,2,2,2))
as.numeric(Most_redundant_order)
j= 1
for (i in rev(Most_redundant_order)) {
  Plot_traits_vs_attributes_highlight(i,stc_highlights[j])
  j = j+1
}

# module categories something
module_categories <- Get_module_categories(kegg_categories_script)
a <- getMetricDistModule(metric = 'PC', module_categories)

# Model module
Module_Names <- RNAseq.data$features$trait_presence_absence[,'3300026282-bin.4'] %>% which(. == TRUE) %>% names

model_bin = '3300026282-bin.4'
margins =c(5,5)
Model_Module_List_All <- Model_Module(RNAseq.data, trait.attributes, model_bin, Module_Names, bkgd.traits)
Plot_Model_Module(Model_Module_List_All, model_bin, Module_Names,margins, sortbygenome='3300026282-bin.4')

# looking at specific bins
rownames(Model_Module_List_All[["Model_Comparison_Matrix"]])
sort(Model_Module_List_All[["Model_Comparison_Matrix"]][15,])

# Write out normalized table to explore downstream
normalized <- RNAseq.data[["table"]]
write_delim(normalized, "~/Desktop/normalized-tbasco-table.tsv", delim="\t")
