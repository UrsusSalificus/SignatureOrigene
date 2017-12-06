#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(ggplot2)

### Boxplot comparison of mean distance to center
# Will input the arguments:
# 1. Percentages obtained through 
# 2. Output file
args <- commandArgs(trailingOnly=TRUE)

output = args[2]

percentages <- read.delim(args[1], header = FALSE)
colnames(percentages) <- c('species', 'factor', 'perc')
dict <- readRDS('../input/dictionary.RData')

levels(percentages$species) <- unname(sapply(levels(as.factor(percentages$species)), function(each) {dict$true[dict$abbrev == each]}))

# We will also prepare the colours
all_colours <- c("#3478f3", "#f7262b", "#10b9cf", "#dab000", "#c841ca", "#00982a", "#f75d21")
# We want the whole genome to be black, so we pick the number of non-whole colour, then add black at the end
colours <- c(all_colours[1:(nlevels(as.factor(percentages$factor)))])

plot_title <- 'Percentages of nucleotides within each factor in the whole genome'


# Invisible points -> enable us to control the color legend to have legend we want
# Note the shape 15 -> square for square in legend...
png(output, width=1600, height=900, units="px")
ggplot(percentages, aes(x = species, y = perc, fill = factor, label = perc)) +
  geom_bar(stat = "identity") +
  ylim(0, 100) +
  labs(title=plot_title, x ="Species", y = "Percentage [%]") +
  theme(
    plot.title = element_text(size = 20, face="bold"),
    axis.title.x = element_text(size = 20),
    axis.text.x  = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    axis.text.y  = element_text(size = 18),
    legend.title = element_text(size=19, face = 'bold'),
    legend.text = element_text(size=18),
    legend.key = element_blank()
  ) +
  scale_fill_manual(name = 'Factors', values = colours, 
                          labels = c(levels(as.factor(percentages$factor))))
dev.off() 


