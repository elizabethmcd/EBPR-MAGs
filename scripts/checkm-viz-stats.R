# Visualiation of completion/contamination of genome bins

library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(plotly)
library(data.table)
library(gtable)

# Read in concatenated lineage text file
lineages = read.table("results/stats/refined-bin-stats.txt", header=TRUE)


# Set significant figures for completion and contamination 
lindf <- as.data.frame(lineages)
lindf$Completeness <- formatC(lindf$completeness, digits=2, format="fg")
lindf$Redundancy <- formatC(lindf$redundancy, digits=2, format="fg")


# Plot
p1 <- lineages %>% ggplot(aes(x=completeness, y=redundancy, color=bin)) + geom_point(size=3)
p2 <- p1 + theme_bw() + theme(legend.position="none") + scale_y_log10()
p2
ggsave("figs/bins-stats-plot.png", p2, width=20, height=10, units=c("cm"))
ggplotly(p2)

# Combine with classification info
ebpr_maned = read.csv("results/mapping/ebpr-rank-abundance-classifications.csv", header=TRUE)
linTax = left_join(lindf, ebpr_maned)

# plot with highest classification colors
fin <- linTax %>% ggplot(aes(x=completeness, y=redundancy, color=Highest_Classf)) + geom_point(size=4) + scale_color_brewer(palette="Paired") + theme_classic()

# save plot
ggsave("figs/bin-stats-colorCoord-plot.png", fin, width=20, height=10, units=c("cm"))
