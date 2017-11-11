#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(ggplot2)

### MultiDimensional Scaling (MDS) analysis of a distance matrix
# Will input the arguments:
# 1. path to the distance matrix

args <- commandArgs(trailingOnly=TRUE)

distance_matrix <- readRDS(args[1])
distance_matrix <- as.matrix(distance_matrix)

# Mean_distance
mean_dist <- rowMeans(distance_matrix)

sorted <- order(mean_dist)
# We will use the 10% first (and thus more similar) windows
cat(sorted[1:round(length(sorted) * 0.1)])
