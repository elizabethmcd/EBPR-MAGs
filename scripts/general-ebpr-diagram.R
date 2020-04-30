library(tidyverse)

ebpr <- read.csv("metadata/general-ebpr-data.csv")

cycle <- ebpr %>% ggplot(aes(x=Minutes)) + geom_line(aes(y=P), colour="#384B60", size=2.5) + geom_line(aes(y=Acetate), colour="#FF6E3C", size=2.5) + geom_line(aes(y=PHA), colour="#205C40", size=2.5) + scale_x_continuous(breaks=seq(0,280, by=50), expand=c(0,0)) + scale_y_continuous(limits=c(0,80), breaks = scales::breaks_pretty(n=10), expand=c(0,0)) + geom_vline(xintercept=100, linetype="dotted", size=1) + theme_classic()

transp <- cycle +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

ggsave(plot=transp, file="~/Desktop/ebpr-cycle-diagram.png", width=8, height=5, units=c("in"), bg="transparent")
