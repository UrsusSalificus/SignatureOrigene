#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(hyperSpec)
library(RColorBrewer)
library(ggplot2)


distance_matrix <- readRDS('/home/titouan/PycharmProjects/Master/Main/FCGR/files/distances/pearson/15000_7/hsap_sample_dist_matrix.RData')
matr <- as.matrix(distance_matrix)
mean_dist <- rowMeans(matr)

fit <- readRDS('/home/titouan/PycharmProjects/Master/Main/FCGR/files/distances/pearson/15000_7/hsap_sample_fit.RData')

cluster <- hclust(distance_matrix, method = "ward.D2")
plot(cluster)

clusterCut <- cutree(cluster, 4)

table(clusterCut)

except <- which(clusterCut == '4')

plot(except)

dat <- data.frame(MDS_1 = fit$points[,1], MDS_2 = fit$points[,2], cluster=clusterCut)

plot(dat$MDS_1, dat$MDS_2)
points(dat$MDS_1[except], dat$MDS_2[except], col = 'red', pch = 16)

plot(1:8846, rep(1, time = length(cluster$height)))
points(except, rep(1, time = length(except)), col = 'red', pch = 16, cex = 1.2)

plot(mean_dist)

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





