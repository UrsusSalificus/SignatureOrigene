#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(data.table)
library(hyperSpec)
library(RColorBrewer)
library(ggplot2)

### Infer recombination rates, using as input the coordinates files produced by the homonyme Python script.
# Will input the arguments:
# 1. path to the input file = coordinates of each of the windows of the species' genome
# 2. path to the RR spline functions of the species
# 3. path to the output table of recombination rate per windows

distance_matrix <- readRDS('/home/titouan/PycharmProjects/Master/Main/FCGR/files/distances/pearson/15000_7/e_coli_dist_matrix.RData')

fit <- cmdscale(as.matrix(distance_matrix),eig=TRUE, k=2) # k is the number of dim

cluster <- hclust(distance_matrix, method = "ward.D2")
plot(cluster)

clusterCut <- cutree(cluster, 5)

data <- data.frame(MDS_1 = fit$points[,1], MDS_2 = fit$points[,2], cluster=clusterCut)

ggplot(data, aes(x = as.numeric(MDS_1), y = as.numeric(MDS_2))) + 
  geom_point(aes( colour = factor(cluster)),size = 2) + 
  labs(title="plot_title", x ="Coordinate 1", y = "Coordinate 2") +
  theme(
    plot.title = element_text(size = 15, face="bold"),
    axis.title.x = element_text(size = 15),
    axis.text.x  = element_text(size = 12),
    axis.title.y = element_text(size = 15),
    axis.text.y  = element_text(size = 12)
  ) +
  scale_colour_gradientn(name = paste("% of", "CDS"), colours = colorRampPalette(rev(brewer.pal(10, "RdYlBu")))(100))


