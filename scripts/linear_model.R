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
distance_matrix <- as.matrix(distance_matrix)

# Mean_distance
mean_dist <- rowMeans(distance_matrix)


# All the features:
RR <- read.csv(args[3], sep = '\t', header = FALSE)
LCR <- read.csv(args[4], sep = '\t', header = FALSE)
CDS <- read.csv(args[5], sep = '\t', header = FALSE)

model <- lm(as.numeric(mean_dist) ~ RR$V2 + LCR$V2 + CDS$V2)
par(mfrow=c(2,2))
plot(model)
summary(model)


