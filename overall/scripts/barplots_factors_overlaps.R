#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(ggplot2)
library(UpSetR)

### Barplot comparison of overlaps
# Will input the arguments:
# 1. Percentages obtained through extract_overlap_factors.py script
# 2. Output file
args <- commandArgs(trailingOnly=TRUE)

output = args[2]

percentages <- read.delim(args[1], sep = '\t', header = FALSE)
colnames(percentages) <- c('species', 'factor', 'perc', 'comparison')

dict <- readRDS('../input/dictionary.RData')
levels(percentages$species) <- unname(sapply(levels(as.factor(percentages$species)), 
                                             function(each) {dict$true[dict$abbrev == each]}))

# We will also prepare the colours
all_colours <- c("#3478f3", "#f7262b", "#10b9cf", "#dab000", "#c841ca", "#00982a", "#f75d21")
# We want the whole genome to be black, so we pick the number of non-whole colour, then add black at the end
colours <- c(all_colours[1:(nlevels(as.factor(percentages$factor)))])

plot_title <- 'Overlapping percentages of the various factors'

png(output, width=1600, height=900, units="px")
ggplot(percentages, aes(x = factor, y = perc, fill = comparison, label = perc)) +
  geom_bar(stat = "identity") +
  facet_wrap(~species) +
  labs(title=plot_title, x ="Factors", y = "Cumulative percentages [%]") +
  theme(
    plot.title = element_text(size = 20, face="bold"),
    axis.title.x = element_text(size = 20),
    axis.text.x  = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    axis.text.y  = element_text(size = 18),
    legend.title = element_text(size=19, face = 'bold'),
    legend.text = element_text(size=20),
    strip.text = element_text(size=20)
  ) +
  scale_fill_manual(name = 'Factors', values = colours, 
                    labels = c(levels(as.factor(percentages$factor))))
dev.off() 
