#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(data.table)

### Manhattan distance computation
# Will input the arguments:
# 1. path to the input file = genomic signature file, each row = 1 region, which will be pariwise compared
# 2. path to the output dist object as RData (with its name in the path)

args <- commandArgs(trailingOnly=TRUE)

feature_table <- fread(args[1], sep = "\t", header = FALSE)

distance_matrix <- dist(feature_table[,2:(ncol(feature_table)-1)], method = 'manhattan')
attr(distance_matrix, "Labels") <- feature_table$V1

# Save as R dist object, to be imported again later on 
saveRDS(distance_matrix, file = args[2])
