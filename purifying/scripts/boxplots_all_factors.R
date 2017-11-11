#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(RColorBrewer)
library(ggplot2)

### Boxplots of mean distances of distances amtrices of FCGRs
# Will input the arguments:
# 1. path to the output image (with its name in the path)
# 3. path to the feature file
# 4. Genomic signature type
# 5. feature type
# 6. window size
# 7. IF USING FCGR : k-mer size

gs = basename(args[4])
window_size <- args[6]
if (gs == 'FCGRs'){
  kmer <- args[7]
}

dist_matrix_path <- paste('files/distances/pearson/', window_size, '_', kmer, '/', sep = '')

all_matrices <- list.files(path = 'files/distances/pearson/15000_7/')
# We will extract the factor from the file names
factors <- unname(sapply(all_matrices, function(each_dist) {strsplit(each_dist, split = '_')[[1]][3]}))

all_raw_means<- lapply(all_matrices, function(each_dist) {
  dist <- readRDS(paste(dist_matrix_path, each_dist, sep =''))
  if (length(dist) > 0) {
    unname(rowMeans(as.matrix(dist)))
  }
})

# We must remove empty distance matrices (for example not enough data to perform pariwise comparison)
to_keep <- ! sapply(all_raw_means, is.null)

factors <- factors[to_keep]
all_raw_means <- all_raw_means[to_keep]

# We will then order by window number 
# TODO: find a way to remove and order at the same time...
ordering <- order(sapply(all_raw_means, length), decreasing = TRUE)

factors <- factors[ordering]
all_raw_means <- all_raw_means[ordering]

# We also need the factor as a vector of same length than their respective distance matrix
factors_as_variable <- sapply(1:length(factors), function(each_dist) {
  rep(factors[each_dist], time = length(all_raw_means[[each_dist]]))
})

dat <- data.frame(factor = unlist(factors_as_variable), mean_dist = unlist(all_raw_means))
boxplot(dat$mean_dist ~ dat$factor)



