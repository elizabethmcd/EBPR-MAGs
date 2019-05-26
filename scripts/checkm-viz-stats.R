# Visualiation of completion/contamination of genome bins

library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(plotly)
library(data.table)
library(gtable)

# Read in concatenated lineage text file
lineages = read.table("Accum-Div/Data/bins/bins-stats.txt")
names(lineages) <- c("Bin", "Classification", "Genome_Length", "Completeness", "Redundancy")


# Set significant figures for completion and contamination 
lindf <- as.data.frame(lineages)
lindf$Completeness <- formatC(lindf$Completeness, digits=2, format="fg")
lindf$Redundancy <- formatC(lindf$Redundancy, digits=2, format="fg")

# Supplement with Phylodist classifications and match bin name
phylodist = read.table("Accum-Div/Data/bins/bins-phylodist-classfs.txt", sep= "\t", header=FALSE)
as.data.frame(phylodist)
names(phylodist) <- c("Bin", "Classification")
cleanphylo = read.table("Accum-Div/Data/bins/bin-cleaned-classifications.txt", sep="\t", header=TRUE)

# Master with full classification and stats
master = left_join(lineages, cleanphylo, by="Bin")
master.ordered = master %>% setorder(-Completeness)
write.table(master.ordered, file="bins-stats-classifications.txt", col.names=TRUE, sep="\t", row.names=FALSE)
complete = master %>% filter(Completeness > 90) 
ordered = complete %>% setorder(-Completeness)
write.table(ordered, file="top-bins-stats-classifications.txt", col.names=TRUE, sep="\t", row.names=FALSE)

# Plot
p1 <- lineages %>% ggplot(aes(x=Completeness, y=Redundancy, color=Bin)) + geom_point(size=2)
p2 <- p1 + theme_bw() + theme(legend.position="none") + scale_y_log10()
p2
ggsave("bins-stats-plot.png", p2, width=20, height=10, units=c("cm"))
ggplotly(p2)

# Above 90% complete chart
ordered <- as.data.frame(ordered)
ordered$Completeness <- format(ordered$Completeness, digits=2, format="fg")
ordered$Redundancy <- format(ordered$Redundancy, digits=1, format="fg")
ordered$Mbp <- ordered$Genome_Length/1000000
ordered$Mbp <- format(ordered$Mbp, digits=3, foramt='fg')
orderfin <- ordered %>% select(-Genome_Length, -Classification.x, -Bin)
colnames(orderfin)[3] <- "Classification"
setcolorder(orderfin, c("Classification", "Completeness", "Redundancy", "Mbp"))
grid.table(orderfin)
g1 <- tableGrob(orderfin, rows=NULL)
g1 <- gtable_add_grob(g1,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                      t = 2, b = nrow(g1), l = 1, r = ncol(g1))
g1 <- gtable_add_grob(g1,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                      t = 1, l = 1, r = ncol(g1))
grid.draw(g1)
final90 <- grid.arrange(g1, nrow=2, as.table=TRUE, heights=c(3,1))
ggsave(file="bins90.png", final90, width=40, height=30, units=c("cm"))

# top 10 chart
top = read.csv("Data/bins/top-10-bins-reordered.csv", sep=",", header=TRUE)
topdf <- as.data.frame(top)
topdf$Completeness <- format(topdf$Completeness, digits=2, format="fg")
topdf$Redundancy <- format(topdf$Redundancy, digits=1, format="fg")
topdf$Mbp <- topdf$Genome_Length/1000000
topdf$Mbp <- format(topdf$Mbp, digits=3, format="fg")
topfin <- topdf %>% select(-Genome_Length)
grid.table(topfin)
g <- tableGrob(topfin, rows=NULL)
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(g))
grid.draw(g)
final <- grid.arrange(g, nrow=2, as.table=TRUE, heights=c(3,1))
ggsave(file="top-bins-stats-table.png", final, width=30, height=30, units=c("cm"))
