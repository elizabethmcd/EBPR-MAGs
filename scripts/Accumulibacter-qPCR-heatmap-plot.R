# Accumulibacter Clade qPCR Abundance heatmap Visualization

library(ggplot2)
library(dplyr)
library(tidyr)
library(gplots)
library(reshape2)
library(viridis)

# qPCR data
data <- read.table("Accum-Flank/Data/Accum-qPCR-Timeseries.csv", sep=',', header=TRUE)
df <- as.data.frame(data)
df$Date <- as.Date(df$Date, "%m/%d/%y")

# Subset by dates for 16S
t1 <- as.Date("2009-01-02")
t2 <- as.Date("2011-12-02")
sub <- df %>% filter(Date>=t1 & Date<=t2)
sub.m <- melt(sub,id.vars="Date",measure.vars=c("Clade.IA", "Clade.IIA"))
sub.m$Date <- as.character(sub.m$Date)

# Full timeseries
sub.f <- melt(df, id.vars="Date", measure.vars=c("Clade.IA", "Clade.IIA"))
sub.f$Date <- as.character(sub.f$Date)

# subset plots 

p1 <- ggplot(sub.m, aes(x=Date, y=variable, fill=value)) + geom_tile(color="white") + scale_x_discrete(labels=sub.m$Date) + scale_fill_gradient(low="lightblue", high="midnightblue")
p2 <- p1 + theme(axis.text.x= element_text(angle=75, hjust=1)) + guides(fill = guide_colorbar(nbin = 4))
p3 <- p2 + guides(fill=FALSE)
ggsave("Accumulibacter-qPCR-heatmap.png", p3, width=20, height=5, units="cm")

# full plots
p4 <- ggplot(sub.f, aes(x=Date, y=variable, fill=value)) + geom_tile(color="white") + scale_x_discrete(labels=sub.f$Date) + scale_fill_gradient(low="lightblue", high="midnightblue")
p5 <- p4 + theme(axis.text.x= element_text(angle=75, hjust=1)) + guides(fill = guide_colorbar(nbin = 4)) 
p6 <- p5 + guides(fill=FALSE)
p6
ggsave("Accumulibacter-qPCR-timeseries.png", p6, width=35, height=5, units="cm")


