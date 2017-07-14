#!/usr/bin/env Rscript

### PCA analysis of a set of FCGR
# Will input the arguments:
# 1. path to the distance matrix file
# 2. path to the feature file
args <- commandArgs(trailingOnly=TRUE)


distance_matrix <- read.table(args[1], sep = "\t", header = FALSE)
# Will remove any empty column
distance_matrix <- distance_matrix[, colSums(is.na(distance_matrix)) == 0]

# Copying the upper diagonal into the lower one
distance_matrix[lower.tri(distance_matrix)] <- t(distance_matrix)[lower.tri(distance_matrix)]

# Finding which column (= which region) as the smallest distance to the other column/region.
# This column/region will then represent the median region.
min_column = which.min(colSums(distance_matrix))
median_region = distance_matrix[-min_column, min_column]


feature <- read.table(args[2], sep = "\t", header = FALSE)
# We take out the median region's feature
feature <- feature[-min_column,]

plot(median_region, feature[,2])

