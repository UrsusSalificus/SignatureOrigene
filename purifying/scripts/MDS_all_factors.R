#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(RColorBrewer)
library(ggplot2)

### MultiDimensional Scaling (MDS) analysis of a distance matrix
# Will input the arguments:
# 1. path to the output image (with its name in the path)
# 2. path to the fit of the distance matrix
# 3. window size
# 4. k-mer size
args <- commandArgs(trailingOnly=TRUE)

output = args[1]
window_size <- args[3]
kmer <- args[4]

# MDS
fit <- readRDS(args[2])

# We now have to fetch from which factor it is from
all_names <- row.names(fit$points)
factors <- unname(sapply(all_names, function(each_line) {
  if (length(strsplit(each_line, '_')[[1]]) == 1){
    # First case, when split we only get one result (no '_') = concatenated factor
    return(each_line)
  } else if (length(strsplit(each_line, '_')[[1]]) == 2) {
    # Second case, when split we get 2 result = whole genome (only the record name)
    return('whole')
  } else {
    # Third case, when split we get 3 result = pure non-concatenated
    return(strsplit(each_line, '_')[[1]][1])
  }
}))

# We will also prepare the colours
all_colours <- c("#53b125", "#0057d2", "#d99b24", "#d940ad", "#cf290e")
# We want the whole genome to be black, so we pick the number of non-whole colour, then add black at the end
colours <- c(all_colours[1:(nlevels(as.factor(factors))-1)], 'black')

factors <- unname(sapply(all_names, function(each_line) {
  if (length(strsplit(each_line, '_')[[1]]) == 1){
    # First case, when split we only get one result (no '_') = concatenated factor
    return(each_line)
  } else if (length(strsplit(each_line, '_')[[1]]) == 2) {
    # Second case, when split we get 2 result = whole genome (only the record name)
    return('whole')
  } else {
    # Third case, when split we get 3 result = pure non-concatenated
    return(strsplit(each_line, '_')[[1]][1])
  }
}))

# We will reduce visibility of the whole genome points:
reduced_vis <- unname(sapply(factors, function (each_line) {
  if (each_line == 'whole'){
    return(TRUE)
  } else {
    return(FALSE)
  }
}))

data <- data.frame(MDS_1 = fit$points[,1], MDS_2 = fit$points[,2], factors = factors, reduced_vis = reduced_vis)

plot_title <- paste('Distance matrices of either pure factor sequences or whole genome,\nwith kmer = ', 
                    kmer, ' and ', window_size, ' bp windows ', sep = '')

png(output, width=900, height=650, units="px")
ggplot(data, aes(x = MDS_1, y = MDS_2, colour = as.factor(factors), alpha = reduced_vis)) + 
  geom_point(size = 2) +
  labs(title=plot_title, x ="Coordinate 1", y = "Coordinate 2") +
  theme(
    plot.title = element_text(size = 20, face="bold"),
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 15),
    legend.title = element_text(size=18),
    legend.text = element_text(size = 15)
  )+
  scale_colour_manual(name = 'Factors', values = colours, 
                      labels = c(levels(as.factor(factors))[-nlevels(as.factor(factors))], 'Whole\ngenome')) +
  scale_alpha_discrete(range = c(1, 0.15), guide=FALSE)
dev.off() 




