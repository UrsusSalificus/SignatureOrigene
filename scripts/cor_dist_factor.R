#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

### Correlation analysis of a a feature in function of distance to median region
# Will input the arguments:
# 1. path to the output image (with its name in the path)
# 2. path to the distance matrix
# 3. path to the feature file
# 4. Genomic signature type
# 5. feature type
# 6. window size
# 7. IF USING FCGR : k-mer size
args <- commandArgs(trailingOnly=TRUE)

output <- args[1]
gs <- basename(args[4])
feature_type <- args[5]
window_size <- args[6]
if (gs == 'FCGRs'){
  kmer <- args[7]
} 

distance_matrix <- readRDS(args[2])
distance_matrix <- as.matrix(distance_matrix)

# Mean_distance
mean_dist <- rowMeans(distance_matrix)

feature <- read.table(args[3], sep = "\t", header = FALSE)

#  Correlation test between distance to median and feature; and parsing the resulting p-value
pv <- round(cor.test(mean_dist, feature[,2], method = 'spearman')$p.value, 6)

if (gs == 'FCGRs'){
  plot_title <- paste('% of', feature_type, 'in function of distance to the median region \n',
                      paste('(FCGRs, k= ', kmer, ', ', window_size, ' bp windows)', '\n', sep = ''),
                      'Correlation test p-value :', pv)

} else {
  plot_title <- paste('% of', feature_type, 'in function of distance to the median region \n',
                      paste('(DFTs, ', window_size, ' bp windows)', '\n', sep = ''),
                      'Correlation test p-value :', pv)
}


png(output, width=900, height=650, units="px")
plot(mean_dist, feature[,2], main=plot_title, xlab='Distance to median region',
     ylab = paste('% of', feature_type), pch = 19)
dev.off() 


