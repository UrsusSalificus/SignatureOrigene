#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(wesanderson)
library(ggplot2)

### MultiDimensional Scaling (MDS) analysis of a distance matrix
# Will input the arguments:
# 1. path to the output image (with its name in the path)
# 2. path to the fit of the distance matrix
# 3. path to the feature file
# 4. Genomic signature type
# 5. factor type
# 6. window size
# 7. sample size
# 8. species
# 9. IF USING FCGR : k-mer size
args <- commandArgs(trailingOnly=TRUE)

output = args[1]
gs = basename(args[4])
factor_type <- args[5]
window_size <- args[6]
sample_size <- args[7]
species <- args[8]
if (gs == 'FCGRs'){
  kmer <- args[9]
} 

# We will need a dictionnary of species abbreviation VS nice species name
dict <- readRDS('../input/dictionary.RData')

# MDS
fit <- readRDS(args[2])
data <- data.frame(MDS_1 = fit$points[,1], MDS_2 = fit$points[,2])

feature = read.csv(args[3], sep = '\t', header = FALSE)

# We will check how much we can round the result and still not lose significant differences:
data$size <- feature[,2]
good_decimal <- 1
round_up <- TRUE
while (round(min(data$size), good_decimal) == 0) {
  good_decimal = good_decimal + 1
  # If the minimum region contain less than 0.1% CDS, stop here
  if (good_decimal == 3){
    break
  }
}
data$size <- round(data$size, good_decimal)

# Plot title: varying in between DFTs and FCGRs
if (gs == 'FCGRs'){
  plot_title <- paste(dict$true[dict$abbrev == species], ': distance matrix of ', sample_size, 
                      ' samples k-mer frequencies,\nwith k-mer = ', 
                      kmer, ' and ', window_size, ' base pairs per sample', sep = '')
} else {
  plot_title <- paste(dict$true[dict$abbrev == species], ': distance matrix of ', sample_size, 
                      ' samples CGR power spectrums,\nwith k-mer = ', 
                      kmer, ' and ', window_size, ' base pairs per sample', sep = '')
}

# Legend: vayring between RR and other factors
if (factor_type == 'RR'){
  colour_legend <- 'Recombination\nrate'
} else {
  colour_legend <- paste("% of", factor_type)
}

png(output, width=1000, height=700, units="px")
ggplot(data, aes(x = MDS_1, y = MDS_2, colour = size)) + 
  geom_point(size = 3.5) + 
  labs(title=plot_title, x ="Coordinate 1", y = "Coordinate 2") +
  theme(
    plot.title = element_text(size = 20, face="bold"),
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 15),
    legend.title = element_text(size=18, face = 'bold'),
    legend.text = element_text(size = 15)
  ) + 
  scale_colour_gradient(name = colour_legend, low = '#1252F1', high = '#e11339')
dev.off() 

