library(dplyr)
library(ggplot2)
library(reshape2)
library(gtable)
library(gridExtra)
library(grid)

# Read in master file of coverage statistics for every bin to each timepoint
covg = read.delim("all-bins-coverage.txt", sep="\t", header=TRUE)

# Classification and genome stats 
binstats = read.delim("McMahon-Lab/EBPR-Projects/EBPR-MAGs/results/refined-bin-stats.txt", sep="\t", header=TRUE)
colnames(binstats)[1] = "ref"
covg = covg %>% filter(filename !="filename")
covg$AvgCov = as.numeric(as.character(covg$AvgCov))
ebpr = left_join(covg, binstats, by.y=ref)
covg$ref = gsub("-", "", covg$ref)
ebprdrop = ebpr %>% filter(meta != "2005-06-14-EBPR")


bubble <- ebprdrop %>% ggplot(aes(x=meta, y=ref)) + geom_point(aes(colour=classification, size=AvgCov)) + theme(panel.background = element_rect(fill="white"), panel.grid.major=element_line(size=.5,colour="grey95", linetype="dotted"), panel.grid.minor=element_line(size=.5, colour="grey95", linetype="dotted"), axis.text.x=element_text(angle=60,size=10,vjust=.5), axis.text.y=element_text(size=7), axis.line=element_blank())

refinedbins = read.delim("results/refined-bin-stats.txt", sep="\t", header=TRUE)
refinedbins$size = refinedbins$size / 1000000
refinedbins$completeness = format(refinedbins$completeness, digits=2, format="fg")
refinedbins$redundancy = format(refinedbins$redundancy, digits=1, format="fg")
bintable = refinedbins[,c(1,2,3,4,5,6,9)]
g1 = tableGrob(bintable, rows=NULL)
g1 <- gtable_add_grob(g1,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                      t = 2, b = nrow(g1), l = 1, r = ncol(g1))
g1 <- gtable_add_grob(g1,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                      t = 1, l = 1, r = ncol(g1))
grid.draw(g1)
final= grid.arrange(g1, nrow=2, as.table=TRUE, heights=c(3,1))
ggsave(file="binstats.png", final, width=50, height=60, units=c("cm"))

p1 <- bintable %>% ggplot(aes(x=completeness, y=redundancy, color=classification)) + geom_point(size=6, alpha=0.5)
p2 <- p1 + theme_bw() + theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15))
p2
ggsave(file="binplot.png", p2, width=35, height=15, units=c("cm"))
ggsave(file="relabundbins.png", bubble, width=40, height=20, units=c("cm"))
