#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(RColorBrewer)
library(ggplot2)

### MultiDimensional Scaling (MDS) analysis of a distance matrix
# Will input the arguments:
# 1. path to the output image (with its name in the path)
# 2. path to the distance matrix
# 3. path to the feature file
# 4. Genomic signature type
# 5. feature type
# 6. window size
# 7. IF USING FCGR : k-mer size
args <- commandArgs(trailingOnly=TRUE)

output = args[1]
gs = basename(args[4])
feature_type <- args[5]
window_size <- args[6]
if (gs == 'FCGRs'){
  kmer <- args[7]
} 

distance_matrix <- readRDS(args[2])

# MDS
fit <- cmdscale(distance_matrix,eig=TRUE, k=2) # k is the number of dim
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

if (gs == 'FCGRs'){
  plot_title <- paste('MDS metrics of FCGRs distance matrix,\nwith kmer = ', 
                      kmer, ' and ', window_size, ' bp windows ', sep = '')
} else {
  plot_title <- paste('MDS metrics of DFTs distance matrix,\nwith ',
                      window_size, ' bp windows distance matrix', sep = '')
}

png(output, width=700, height=500, units="px")
ggplot(data, aes(x = MDS_1, y = MDS_2, colour = size)) + 
  geom_point(size = 2) + 
  labs(title=plot_title, x ="Coordinate 1", y = "Coordinate 2") +
  theme(
    plot.title = element_text(size = 15, face="bold"),
    axis.title.x = element_text(size = 15),
    axis.text.x  = element_text(size = 12),
    axis.title.y = element_text(size = 15),
    axis.text.y  = element_text(size = 12)
  ) +
  scale_colour_gradientn(name = paste("% of", feature_type), colours = colorRampPalette(rev(brewer.pal(10, "RdYlBu")))(100))
dev.off() 

