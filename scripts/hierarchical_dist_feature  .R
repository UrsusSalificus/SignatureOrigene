#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(hyperSpec)
library(RColorBrewer)
library(ggplot2)


FCGRs <- read.table("files/FCGRs/5000_4/s_cerevisiae_FCGRs.txt", sep = "\t", header = FALSE)
distance_matrix <- pearson.dist(FCGRs[,c(-1, -ncol(FCGRs))])

# MDS
fit <- cmdscale(distance_matrix,eig=TRUE, k=2) # k is the number of dim

cluster <- hclust(distance_matrix, method = "ward.D2")
plot(cluster)

clusterCut <- cutree(cluster, 5)

data <- data.frame(MDS_1 = fit$points[,1], MDS_2 = fit$points[,2], cluster=clusterCut)

feature = read.csv("../files/features/5000_CDS/s_cerevisiae.txt", sep = '\t', header = FALSE)

# We will check how much we can round the result and still not lose significant differences:
data$size <- feature[,2]

ggplot(data, aes(x = MDS_1, y = MDS_2)) + 
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





