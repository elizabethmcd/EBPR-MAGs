library(TbasCO)

# File of raw counts and annotations
ebpr_counts = file.path("data/raw-ebpr-counts-data.csv")

# Normalization
normalization.features <- list('no_feature'   = c(9159700, 4459877, 9826273, 8171512, 9542765, 10522313),
  'ambiguous'= c(3940698, 2023389, 4675033, 3308789, 6446272, 5966543),
  'library_size' = c(234232896, 183166236, 228746720, 198024002, 231567992, 259156166),
  'not_aligned'  = c(0, 0, 0, 0, 0, 0)
)

# RNA seq object
RNAseq.data <- Pre_process_input(ebpr_counts, annotation.db.path, normalize.method=T, filter.method='MAD', filter.low.coverage = T, normalization.features)

# Find bins that made above cutoff threshhold
as.data.frame(rownames(Model_Module_Matrix[["Model_Sig_Matrix"]]))

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
bkgd.individual <- Individual_Annotation_Background(RNAseq.data,N= 1000, metrics = distance.metrics, threads = 2)

bkgd.individual.Zscores <- Calc_Z_scores(bkgd.individual, distance.metrics)

bkgd.traits <- Random_Trait_Background(RNAseq.data, bkgd.individual.Zscores,N=10000,metrics=distance.metrics,threads=2)

pairwise.distances  <- Calc_Pairwise_Annotation_Distance(RNAseq.data, RNAseq.data$features$annotation.db, distance.metrics,bkgd.individual.Zscores, show.progress = T,threads = 4)

trait.attributes <- Identify_Trait_Attributes(RNAseq.data = RNAseq.data,
                                              pairwise.distances = pairwise.distances,
                                              annotation.db = RNAseq.data$features$annotation.db,threads = 2)

trait.attributes.pruned <- Prune_Trait_Attributes(trait.attributes, bkgd.traits, 
                                                  RNAseq.data, 
                                                  p.threshold = 0.05,
                                                  pairwise.distances = pairwise.distances, 
                                                  bkgd.individual.Zscores = bkgd.individual.Zscores, 
                                                  annotation.db = RNAseq.data$features$annotation.db,
                                                  trait_presence_absence = RNAseq.data$features$trait_presence_absence)

sbs.trait.attributes <- Traitattributes_To_Sbsmatrix(trait.attributes.pruned,RNAseq.data$features$bins)

# plot background of individual genes
Plot_Background_Individual_Genes(bkgd.individual.Zscores)

# comparison of metrics plots
Plot_Metric_Comparison <- function(bkgd.individual) {
  t_test_KO_random_pearson <- t.test(bkgd.individual$`Random Annotated Genes`$PC,
                                     bkgd.individual$`Random Genes`$PC,
                                     alternative = "less"
  ) # x > y (NULL)
  
  t_test_KO_random_euclidean <- t.test(bkgd.individual$`Random Annotated Genes`$NRED,
                                       bkgd.individual$`Random Genes`$NRED,
                                       alternative = "greater"
  ) # x > y (NULL)
  
  par(
    mfrow = c(2, 2),
    mar = c(3, 3, 3, 1)
  )
  # plot 1
  plot(density(bkgd.individual$`Random Genes`$PC,
               adjust = 2,
               na.rm = TRUE
  ),
  ylim = c(0, 1),
  xlab = "",
  ylab = "",
  main = ""
  )
  points(density(bkgd.individual$`Random Annotated Genes`$PC, adjust = 2),
         typ = "l",
         col = "blue"
  )
  mtext(paste("p-value = ", signif(t_test_KO_random_pearson$p.value, 2)),
        side = 3,
        col = "blue",
        padj = 2,
        cex = 0.75
  )
  title(
    ylab = "Density",
    line = 2,
    cex.lab = 1
  )
  title(
    xlab = "PC",
    line = 2,
    cex.lab = 1
  )
  
  # plot 2
  plot(density(bkgd.individual$`Random Annotated Genes`$PC,
               adjust = 2
  ),
  ylim = c(0, 1),
  xlab = "",
  ylab = "",
  main = " "
  )
  points(density(bkgd.individual$`Genes with the same annotation`$PC,
                 adjust = 2
  ),
  typ = "l",
  col = "red"
  )
  mtext(paste("p-value = ", signif(t_test_KO_random_pearson$p.value, 2)),
        side = 3,
        col = "red",
        padj = 2,
        cex = 0.75
  )
  title(
    ylab = "Density",
    line = 2,
    cex.lab = 1
  )
  title(
    xlab = "PC",
    line = 2,
    cex.lab = 1
  )
  
  # plot 3
  plot(density(bkgd.individual$`Random Genes`$NRED,
               adjust = 2
  ),
  typ = "l",
  ylim = c(0, 1.25),
  xlab = "",
  ylab = "",
  main = ""
  )
  points(density(bkgd.individual$`Random Annotated Genes`$NRED, adjust = 2),
         typ = "l",
         col = "blue"
  )
  title(
    ylab = "Density",
    line = 2,
    cex.lab = 1
  )
  title(xlab = "NRED", line = 2, cex.lab = 1)
  mtext(paste("p-value = ", signif(t_test_KO_random_euclidean$p.value, 2)),
        side = 3, col = "blue", padj = 2, cex = .75
  )
  
  # plot 4
  plot(density(bkgd.individual$`Random Annotated Genes`$NRED,
               adjust = 2
  ),
  typ = "l",
  ylim = c(0, 1.25),
  xlab = "",
  ylab = "",
  main = ""
  )
  points(density(bkgd.individual$`Random Annotated Genes`$NRED,
                 adjust = 2
  ),
  typ = "l",
  col = "red"
  )
  title(
    ylab = "Density",
    line = 2,
    cex.lab = 1
  )
  title(
    xlab = "NRED",
    line = 2,
    cex.lab = 1
  )
  title(" \n\nComparison of random & functional \n pairwise comparisons",
        outer = TRUE
  )
  mtext(paste(
    "p-value = ",
    signif(t_test_KO_random_euclidean$p.value, 2)
  ),
  side = 3,
  col = "red",
  padj = 2,
  cex = .75
  )
}
Plot_Metric_Comparison(bkgd.individual)

# Redundant Traits Function
Plot_Redundancy_Traits <- function(RNAseq.data) {
  library(ggplot2)
  library(magrittr)
  
  ta.pa <- apply(RNAseq.data$features$trait_presence_absence, 1, function(row) {
    return(sum(row))
  })
  
  ta.pa <- ta.pa[which(ta.pa != 0)]
  sort(ta.pa) %>% barplot(xaxt='n', ylim = c(0,25))
  
}
Plot_Redundancy_Traits(RNAseq.data)

# Background of individual genes
Plot_Background_Modules(bkgd.traits)

# Model module
Module_Names <- RNAseq.data$features$trait_presence_absence[,'3300009517-bin.1'] %>% which(. == T) %>% names
margins =c(5,23)
model_bin = '3300026284-bin.9'
Model_Module_List_All <- Model_Module(RNAseq.data, trait.attributes, model_bin, bkgd.traits)
Plot_Model_Module(Model_Module_List_All, model_bin, Module_Names,margins, sortbygenome='3300009517-bin.1')


# looking at specific bins
rownames(Model_Module_Matrix[["Model_Comparison_Matrix"]])
sort(Model_Module_Matrix[["Model_Comparison_Matrix"]][13,])

# Figures of functional redundancy and linear regression
Calc_TnA_redundancy <- function() {
  point.matrix <- matrix(ncol=4, nrow=0)
  
  t.pa <- RNAseq.data$features$trait_presence_absence
  
  combinations <- combn(RNAseq.data$features$bins, 2, simplify = F)
  
  for (pair in combinations){
    try({
      A <- pair[1] %>% as.character
      B <- pair[2] %>% as.character
      
      traits.A <- t.pa[, A] %>% which(. == T) %>% names
      traits.B <- t.pa[, B] %>% which(. == T) %>% names
      overlap.traits <- intersect(traits.A, traits.B) %>% length
      
      attributes.A <- which(sbs.trait.attributes[A,]  == 1) %>% names
      print(B)
      attributes.B <- which(sbs.trait.attributes[B,]  == 1) %>% names
      overlap.attributes <- intersect(attributes.A, attributes.B) %>% length
      
      point.matrix <- rbind(point.matrix,
                            c(A, B,
                              overlap.traits,
                              overlap.attributes))
    })
    
  }
  colnames(point.matrix) <- c("A","B", 'traits', 'attributes')
  
  point.df <- data.frame(point.matrix, stringsAsFactors = F)
  point.df$traits %<>% as.numeric
  point.df$attributes %<>% as.numeric
  all_bins <- unique(c(point.df[,1],point.df[,2]))
  sum_traits <- NULL
  sum_attributes <- NULL
  for (i in as.numeric(all_bins)) {
    bin_rows <- which(point.df[,1] == i | point.df[,2] == i)
    sum_traits <- c(sum_traits, sum(point.df[bin_rows,3]))
    sum_attributes <- c(sum_attributes,sum(point.df[bin_rows,4]))
  }
  return(cbind(all_bins, sum_traits, sum_attributes))
}         

TnA_redundancy <- Calc_TnA_redundancy()
# Sort based on attributes
Most_redundant_order <- TnA_redundancy[order(TnA_redundancy[,3]),1]
# Color by taxonomy
stc<- string.to.colors(genome.taxonomy)
# Color by taxonomy in the order based on redundancy
stc_highlights <- stc[match(rev(as.numeric(Most_redundant_order)),names(stc))]
par(mfrow=c(4,5),mar=c(2,2,2,2))
plot(as.numeric(TnA_redundancy[,3])~as.numeric(TnA_redundancy[,2]), col=stc, xlim=c(500,1100), ylim=c(100,400), pch = 19, xlab = "Attributes", ylab = "Traits")

as.numeric(Most_redundant_order)
j= 1
for (i in rev(as.numeric(Most_redundant_order))) {
  Plot_traits_vs_attributes_highlight(i,stc_highlights[j])
  j = j+1
}

calcmfrow <- function(x) {
  temp <- sqrt(x)
  
  if ((temp %% floor(temp)) == 0) {
    return(c(temp, temp))
  } else if ((temp %% floor(temp)) < 0.5) {
    return(c(floor(temp), ceiling(temp)))
  } else {
    return(c(ceiling(temp), ceiling(temp)))
  }
}


Plot_Trait_Attribute_Expression <- function(trait.attribute,
                                            trait.attributes.pruned,
                                            RNAseq.data) {
  trait.attribute.s <- unlist(strsplit(x = trait.attribute, split = "[.]"))
  
  trait.annotations <- RNAseq.data$features$annotation.db$module.dict[[trait.attribute.s[1]]]
  attribute.genomes <- trait.attributes[[trait.attribute.s[1]]][[trait.attribute.s[2]]]$genomes
  
  par(mfrow = calcmfrow(length(trait.annotations)))
  
  for (annotation in trait.annotations) {
    annotation.expression <- RNAseq.data$table[which(RNAseq.data$table$Annotation == annotation &
                                                       RNAseq.data$table$Bin %in% attribute.genomes), ]
    if (nrow(annotation.expression) == 1) {
      expression <- annotation.expression[RNAseq.data$features$rank.columns]
      plot(as.character(expression),
           type = "l", ylab = "Expression value", xlab = "Points in time",
           col = "black", lwd = "3", ylim = c(0, max(expression)), main = annotation
      )
      next()
    }
    mean.cols <- apply(
      annotation.expression[RNAseq.data$features$rank.columns],
      2, mean
    )
    sd.cols <- apply(
      annotation.expression[RNAseq.data$features$rank.columns],
      2, sd
    ) / 2
    
    mean.psd <- mean.cols + sd.cols
    mean.msd <- mean.cols - sd.cols
    
    max.val <- max(mean.psd)
    if (!is.nan(max.val)) {
      plot(mean.cols,
           type = "l", ylab = "Expression value", xlab = "Points in time",
           col = "black", lwd = "3", ylim = c(0, 1), main = annotation
      )
      points(mean.psd, type = "l", col = "blue", lty = 2)
      points(mean.msd, type = "l", col = "red", lty = 2)
    } else {
      plot(c(0, 20, 60, 100, 100, 100),
           type = "n", main = annotation, ylab = "Expression value",
           xlab = "Points in time"
      )
    }
  }
}
Plot_Trait_Attribute_Expression()
Plot_Trait_Attribute_Expression(trait.attribute='M00009_756', trait.attributes$trait.attribute, RNAseq.data)
Plot_Trait_Attribute_Expression()