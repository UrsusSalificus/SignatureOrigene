#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(ggplot2)

### Boxplot comparison of mean distance to center
# Will input the arguments:
# 1. path to the output image (with its name in the path)
# 2. window size
# 3. k-mer size
# 4. species 
args <- commandArgs(trailingOnly=TRUE)

output = args[1]
window_size <- args[2]
kmer <- args[3]
species <- args[4]

distance_directory <- paste("files/distances/pearson/", window_size, '_', kmer, '/', sep ='')

distance_matrices <- Sys.glob(paste(distance_directory, species, '*', sep = ''))
distance_matrices <- lapply(distance_matrices, function(each_matrix) {
  distance_matrix <- readRDS(each_matrix)
  rowMeans(as.matrix(distance_matrix))
})
# TODO -> have it ot be depandant on the file (ex look for whole to treat it differently)
# Then have it to take out the redunduncy (you only need the column of the factor, then only the row of the factor !)
# For whole genome -> can find the center agian by doing :
sorted <- order(mean_dist)
# We will use the 10% first (and thus more similar) windows
cat(sorted[1:round(length(sorted) * 0.1)])



