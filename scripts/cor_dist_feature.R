#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

### Correlation analysis of a a feature in function of distance to median region
# Will input the arguments:
# 1. path to the output image (with its name in the path)
# 2. NOT USED IN THIS SCRIPT
# 3. path to the distance matrix
# 4. path to the feature file
# 5. feature type
# 6. window size
# 7. Genomic signature type
# 8. IF USING FCGR : k-mer size
args <- commandArgs(trailingOnly=TRUE)

output <- args[1]
gs <- basename(args[5])
feature_type <- args[6]
window_size <- args[7]
if (gs == 'FCGRs'){
  kmer <- args[8]
} 

distance_matrix <- read.table(args[3], sep = "\t", header = FALSE)
# Will remove any empty column
distance_matrix <- distance_matrix[, colSums(is.na(distance_matrix)) == 0]

# Copying the upper diagonal into the lower one
distance_matrix[lower.tri(distance_matrix)] <- t(distance_matrix)[lower.tri(distance_matrix)]

# Finding which column (= which region) as the smallest distance to the other column/region.
# This column/region will then represent the median region.
min_column = which.min(colSums(distance_matrix))
median_region = distance_matrix[-min_column, min_column]


feature <- read.table(args[4], sep = "\t", header = FALSE)
# We take out the median region's feature
feature <- feature[-min_column,]

# Quick correlation test and parsing the resulting p-value
pv <- round(cor.test(median_region, feature[,2])$p.value, 6)

if (gs == 'FCGRs'){
  plot_title <- paste('% of', feature_type, 'in function of distance to the median region \n',
                      paste('(FCGRs, k= ', kmer, ', ', window_size, ' bp windows)', '\n', sep = ''),
                      'Correlation test p-value :', pv)

} else {
  plot_title <- paste('% of', feature_type, 'in function of distance to the median region \n',
                      paste('(DFTs, ', window_size, ' bp windows)', '\n', sep = ''),
                      'Correlation test p-value :', pv)
}


png(output, width=700, height=500, units="px")
plot(median_region, feature[,2], main=plot_title, xlab='Distance to median region',
     ylab = paste('% of', feature_type), pch = 19)
dev.off() 



