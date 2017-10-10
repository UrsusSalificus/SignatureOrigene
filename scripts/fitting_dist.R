#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

### Fitting of distance matrix
# Will input the arguments:
# 1. path to the distance matrix
# 2. path to the output fit as RData (with its name in the path)

args <- commandArgs(trailingOnly=TRUE)

distance_matrix <- readRDS(args[1])

# MDS
fit <- cmdscale(distance_matrix,eig=TRUE, k=2) # k is the number of dim

saveRDS(fit, file = args[2])