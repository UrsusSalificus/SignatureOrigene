#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

### Correlation analysis of a a feature in function of distance to median region
# Will input the arguments:
# 1. path to the distance matrix file
# 2. path to the feature file
# 3. feature as displayed in the title and y axis
# 4. distance matrix information (type of genomic signature/type of distance/window size)
# 5. path to the output image (with its name in the path)
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

# Quick correlation test and parsing the resulting p-value
pv <- round(cor.test(median_region, feature[,2])$p.value, 6)

plot_title <- paste(args[3], 'in function of distance to the median region \n', 
                    paste('(', args[4], ')', '\n', sep = ''), 
                    'Correlation test p-value :', pv)

png(args[5], width=700, height=500, units="px")
plot(median_region, feature[,2], main=plot_title, xlab='Distance to median region', ylab = args[3], pch = 19)
dev.off() 



