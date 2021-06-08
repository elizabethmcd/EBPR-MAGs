# R1R2 MAGs Functional Summaries
library(tidyverse)
library(reshape2)
library(grid)
library(gridExtra)

summaries <- read.csv("results/2013_binning/annotations/summaries/cleaned-matrix.csv")
metadata <- read.csv("results/R1R2-EBPR-MAGs-table.csv") %>% 
  select(Bin, Code, Classification)

colnames(summaries)[1] <- c("Bin")
metadata$Bin <- gsub("-", "_", metadata$Bin)
summary_table <- left_join(metadata, summaries)

write.csv(summary_table, "results/2013_binning/annotations/summaries/metabolic_summaries_table.csv", quote=FALSE, row.names=FALSE)

# selected traits (P,N,PHA,S)
selectedTraits <- read.csv("results/2013_binning/annotations/summaries/metabolic_summaries_table_select_traits.csv")
colnames(selectedTraits) <- c("Bin", "Code", "Classification", "hao", "pmoA", "pmoB", "pmoC", "nirB", "nirD", "nirK", "narGZ", "nosZ", "nosD", "aprA", "aprB", "sat", "dsrA", "dsrB", "nifD", "nifH", "nifK", "phbB", "phbA", "phbC", "phaZ", "ppk1", "pstB", "pstC", "pstA", "phoU", "pstS", "pit", "soxA", "soxX", "soxB", "soxC", "soxY","soxZ", "soxD")

traits <- selectedTraits %>% 
  select(-pmoA, -pmoB, -pmoC, -aprA, -aprB, -sat, -dsrA, -dsrB)

traits$pstABCS <- (traits$pstA + traits$pstB + traits$pstC + traits$pstS) / 4
data.frame(apply(traits, 2, function(x) ifelse(x<1,0,x)))

write.csv(traits, "results/2013_binning/annotations/summaries/trait_summaries.csv", quote=FALSE, row.names = FALSE)

traits_modf <- read.csv("results/2013_binning/annotations/summaries/trait_table.csv")

# split trait groups
nitrogen  <- traits_modf %>% 
  select(Code, hao, nirB, nirD, nirK, narGZ, nosZ, nosD, nifD, nifH, nifK)
sulfur <- traits_modf %>% 
  select(Code, soxA, soxX, soxB, soxC, soxY, soxZ, soxD)
phosphorus <- traits_modf %>% 
  select(Code, ppk1, pstABCS, phoU, pit)
pha <- traits_modf %>% 
  select(Code, phbB, phbA, phbC, phaZ)

bin_order <- c('AUS1',
               'PHYC1',
               'PHYC2',
               'TET1',
               'TET2',
               'LEU1',
               'LEU2',
               'SAL1',
               'NANO1',
               'PROP1',
               'PROP3',
               'PROP2',
               'FIMBRI1',
               'BAC1',
               'BAC2',
               'CHIT1',
               'CHIT2',
               'SAP1',
               'SAP2',
               'LEAD1',
               'RUN1',
               'FLAVO1',
               'CHRYS1',
               'BAC3',
               'IGNAVI1',
               'RTHERM1',
               'ANAER1',
               'HERP1',
               'OBS1',
               'FUSI1',
               'GEMMA1',
               'SACCH1',
               'ALPHA1',
               'CAED1',
               'BREV1',
               'CAULO1',
               'HYPHO1',
               'REYR1',
               'REYR2',
               'ANDERS1',
               'BEIJ1',
               'BEIJ2',
               'BEIJ4',
               'BEIJ3',
               'PHREA1',
               'RHIZO1',
               'RHIZO2',
               'RHIZO3',
               'RHODO1',
               'RHODO2',
               'RHODO3',
               'RICK1',
               'SPHING1',
               'ALIC1',
               'OTTO1',
               'OTTO2',
               'RAM1',
               'RUBRI1',
               'VITREO1',
               'CAPIA',
               'CAPIIA',
               'ZOO1',
               'LEG1',
               'LUTEI1',
               'PSEUDO1',
               'PSEUDO2')

# All traits
traits_table <- traits_modf %>% 
  select(Code, hao, nirB, nirD, nirK, narGZ, nosZ, nosD, nifD, nifH, nifK, pstABCS, phoU, ppk1, pit, phbA, phbB, phbC, phaZ, soxA, soxA, soxB, soxC, soxD, soxX, soxY, soxZ)

traits.melted <- melt(traits_table, id.vars="Code") %>% mutate(Code=factor(Code), Code=factor(Code, levels=c(bin_order)))

all_traits_plot <- traits.melted %>% ggplot(aes(x=Code, y=fct_rev(variable), fill=value)) + geom_tile(color="black", size=0.5, aes(width=1, height=1)) + coord_fixed() + theme(axis.text.x.bottom = element_text(angle=80, hjust=0.95, face="bold"), axis.text.y.left = element_text(face="italic"), axis.ticks.y=element_blank())
all_traits_plot

ggsave("figs/all-traits-heatmap-PNS_PHA.pdf", all_traits_plot, width=25, height=10, units=c("cm"))

# Nitrogen
nitrogen.melted <- melt(nitrogen, id.vars="Code") %>% 
  mutate(Code=factor(Code), Code=factor(Code, levels=c(bin_order)))

nitrogen_plot <- nitrogen.melted %>% ggplot(aes(x=fct_rev(Code), y=variable, fill=value)) + geom_tile(color='black', size=0.5, aes(width=1,height=1)) + coord_fixed() + scale_fill_gradient(low="white", high="blue3") + theme(panel.grid=element_blank(), panel.border=element_blank(), plot.margin=unit(c(0,0,0,0), "cm")) + scale_x_discrete(position="top", expand=c(0,0)) + theme(axis.text.x.top=element_text(angle=85, hjust=0, face="italic"), axis.ticks.x=element_blank()) + labs(x=NULL, y=NULL) +  scale_y_discrete(expand=c(0,0))

nitrogen_clean <- nitrogen_plot + theme(legend.position="none")
nitrogen_clean

# Sulfur
sulfur.melted <- melt(sulfur, id.vars="Code") %>% 
  mutate(Code=factor(Code), Code=factor(Code, levels=c(bin_order)))

sulfur_plot <- sulfur.melted %>% ggplot(aes(x=variable, y=fct_rev(Code), fill=value)) + geom_tile(color='black', size=0.5, aes(width=1,height=1)) + coord_fixed() + scale_fill_gradient(low="white", high="darkorchid4") + theme(panel.grid=element_blank(), panel.border=element_blank(), plot.margin=unit(c(0,0,0,0), "cm")) + scale_x_discrete(position="top", expand=c(0,0)) + theme(axis.text.x.top=element_text(angle=85, hjust=0, face="italic"), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) + labs(x=NULL, y=NULL) +  scale_y_discrete(expand=c(0,0))

sulfur_clean <- sulfur_plot + theme(legend.position="none")

# Phosphorus
phosph.melted <- melt(phosphorus, id.vars="Code") %>%
  mutate(Code=factor(Code), Code=factor(Code, levels=c(bin_order)))

phosph_plot <- phosph.melted %>% ggplot(aes(x=variable, y=fct_rev(Code), fill=value)) + geom_tile(color='black', size=0.5, aes(width=1,height=1)) + coord_fixed() + scale_fill_gradient(low="white", high="violetred4") + theme(panel.grid=element_blank(), panel.border=element_blank(), plot.margin=unit(c(0,0,0,0), "cm")) + scale_x_discrete(position="top", expand=c(0,0)) + theme(axis.text.x.top=element_text(angle=85, hjust=0, face="italic"), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) + labs(x=NULL, y=NULL) +  scale_y_discrete(expand=c(0,0)) 

phosph_clean <- phosph_plot + theme(legend.position="none")

# PHA
pha.melted <- melt(pha, id.vars="Code") %>% 
  mutate(Code=factor(Code), Code=factor(Code, levels=c(bin_order)))

pha_plot <- pha.melted %>% ggplot(aes(x=variable, y=fct_rev(Code), fill=value)) + geom_tile(color='black', size=0.5, aes(width=1,height=1)) + coord_fixed() + scale_fill_gradient(low="white", high="tomato2") + theme(panel.grid=element_blank(), panel.border=element_blank(), plot.margin=unit(c(0,0,0,0), "cm")) + scale_x_discrete(position="top", expand=c(0,0)) + theme(axis.text.x.top=element_text(angle=85, hjust=0, face="italic"), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) + labs(x=NULL, y=NULL) +  scale_y_discrete(expand=c(0,0)) 

pha_clean <- pha_plot + theme(legend.position="none")

n.tmp = ggplot_build(nitrogen_clean)
s.tmp = ggplot_build(sulfur_clean)
p.tmp = ggplot_build(phosph_clean)
b.tmp = ggplot_build(pha_clean)
g1 = ggplot_gtable(n.tmp) ; g2=ggplot_gtable(s.tmp) ; g3=ggplot_gtable(p.tmp) ; g4=ggplot_gtable(b.tmp)
n1=length(n.tmp[["panel"]][["ranges"]][[1]][["x.major_source"]])
n2=length(s.tmp[["panel"]][["ranges"]][[1]][["x.major_source"]])
n3=length(p.tmp[["panel"]][["ranges"]][[1]][["x.major_source"]])
n4=length(b.tmp[["panel"]][["ranges"]][[1]][["x.major_source"]])
g=cbind(g1,g2,g3,g4, size="first")
ggsave(file="~/Desktop/test.png", plot=g, width=45, height=26, units=c("cm"))
