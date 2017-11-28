#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(ggplot2)

### Correlation analysis of a a feature in function of distance to median region
# Will input the arguments:
# 1. path to the output image (with its name in the path)
# 2. path to the distance matrix
# 3. path to the feature file
# 4. Genomic signature type
# 5. factor type
# 6. window size
# 7. sample size
# 8. species
# 9. IF USING FCGR : k-mer size
args <- commandArgs(trailingOnly=TRUE)

output <- args[1]
gs <- basename(args[4])
factor_type <- args[5]
window_size <- args[6]
sample_size <- args[7]
species <- args[8]
if (gs == 'FCGRs'){
  kmer <- args[9]
} 

# We will need a dictionnary of species abbreviation VS nice species name
dict <- readRDS('../input/dictionary.RData')

distance_matrix <- readRDS(args[2])
distance_matrix <- as.matrix(distance_matrix)

# Mean_distance
mean_dist <- rowMeans(distance_matrix)

factor <- read.table(args[3], sep = "\t", header = FALSE)
dat <- data.frame(mean_dist = mean_dist, factor=factor[,2])

#  Correlation test between distance to median and factor; and parsing the resulting p-value
pv <- round(cor.test(mean_dist, dat$factor, method = 'spearman')$p.value, 5)
r <- round(cor.test(mean_dist, dat$factor, method = 'spearman')$estimate, 2)

# Legend: vayring between RR and other factors
if (factor_type == 'RR'){
  factor_legend <- 'Recombination rate'
} else {
  factor_legend <- paste("% of", factor_type)
}

# Plot title: varying in between DFTs and FCGRs
if (gs == 'FCGRs'){
  plot_title <- paste(dict$true[dict$abbrev == species], ': ', factor_legend, ' in function of distance to the median region\n',
                      '(', sample_size, ' FCGRs, k= ', kmer, ', ', window_size, 
                      ' bp windows, correlation test p-value: ', pv, ' and rho: ', r, ')', sep = '')
  
} else {
  plot_title <- paste(dict$true[dict$abbrev == species], ':', factor_legend, 'in function of distance to the median region\n',
                      '(', sample_size, ' power spectrums, k= ', kmer, ', ', window_size, 
                      ' bp windows, correlation test p-value: ', pv, ' and rho: ', r, ')', sep = '')
}

png(output, width=1000, height=700, units="px")
ggplot(dat, aes(x = mean_dist, y = factor)) + 
  geom_point(size = 3.5) +
  geom_smooth(method = lm,  linetype = "dashed", color = "darkred", fill = "black") +
  labs(title=plot_title, x = "Distance to median region", y = factor_legend) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 15)
  )
dev.off() 


