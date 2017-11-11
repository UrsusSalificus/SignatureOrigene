#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(ggplot2)
library(dbscan)

### MultiDimensional Scaling (MDS) analysis of a distance matrix
# Will input the arguments:
# 1. path to the output image (with its name in the path)
# 2. path to the distance matrix
# 3. path to the fit of the distance matrix
# 4. window size
# 5. k-mer size
args <- commandArgs(trailingOnly=TRUE)

output = args[1]
window_size <- args[4]
kmer <- args[5]
dict <- readRDS('../input/dictionary.RData')

# Distance matrix
distance_matrix <- readRDS(args[2])
# MDS
fit <- readRDS(args[3])

# We now have to fetch from which factor it is from, and from which species
all_names <- labels(distance_matrix)
species <- unname(sapply(all_names, function(each_line) {
  species_name <- paste(strsplit(each_line, '_')[[1]][1], strsplit(each_line, '_')[[1]][2], sep = '_')
  return(species_name)
}))
factors <- unname(sapply(all_names, function(each_line) {
  factor <- strsplit(each_line, '_')[[1]][3]
  return(factor)
}))

# We will also prepare the colours
all_colours <- c("#53b125", "#0057d2", "#d99b24", "#d940ad", "#cf290e")
# We want the whole genome to be black, so we pick the number of non-whole colour, then add black at the end
colours <- c(all_colours[1:(nlevels(as.factor(factors)))])

data <- data.frame(MDS_1 = fit$points[,1], MDS_2 = fit$points[,2], factors = factors, species = species)

plot_title <- paste('Distance matrices of H.sapiens/M.musculus pure factor sequences,\nwith kmer = ', 
                    kmer, ' and ', window_size, ' bp windows ', sep = '')

species_nice_names <- unname(sapply(levels(as.factor(species)), function (each) {
  return(dict$true[dict$abbrev == each])
}))

png(output, width=1200, height=800, units="px")
ggplot(data, aes(x = MDS_1, y = MDS_2, shape = as.factor(factors), colour = as.factor(species))) + 
  geom_point(size = 2, alpha = 1) +
  labs(title=plot_title, x ="Coordinate 1", y = "Coordinate 2") +
  theme(
    plot.title = element_text(size = 20, face="bold"),
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 15),
    legend.title = element_text(size=18, face = 'bold'),
    legend.text = element_text(size = 15),
    legend.direction = "horizontal",
    legend.position="bottom",
    legend.box = "vertical"
  ) +
  scale_shape_discrete(name = 'Factors', labels = levels(as.factor(factors))) +
  scale_colour_manual(name = 'Species', values = colours, labels = species_nice_names) 
dev.off() 



# TESTING PART
cluster2 <- kmeans(distance_matrix, centers = 2)

table(species, cluster2$cluster)
table(factors, cluster2$cluster)

plot_title <- paste('Clustering of H.sapiens/M.musculus pure factor sequences,\nwith kmer = ', 
                    kmer, 'kmean = 2 and ', window_size, ' bp windows ', sep = '')

png('hsap_sample_mmus_sample_2cluster.png', width=1200, height=800, units="px")
ggplot(data, aes(x = MDS_1, y = MDS_2, shape = as.factor(factors), colour = as.factor(cluster2$cluster))) + 
  geom_point(size = 2, alpha = 1) +
  labs(title=plot_title, x ="Coordinate 1", y = "Coordinate 2") +
  theme(
    plot.title = element_text(size = 20, face="bold"),
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 15),
    legend.title = element_text(size=18, face = 'bold'),
    legend.text = element_text(size = 15),
    legend.direction = "horizontal",
    legend.position="bottom",
    legend.box = "vertical"
  ) +
  scale_shape_discrete(name = 'Factors', labels = levels(as.factor(factors))) +
  scale_colour_manual(name = 'Clustering', values = rev(colours)[4:5]) 
dev.off()

cluster5 <- kmeans(distance_matrix, centers = 5)

table(species, cluster5$cluster)
table(factors, cluster5$cluster)

plot_title <- paste('Clustering of H.sapiens/M.musculus pure factor sequences,\nwith kmer = ', 
                    kmer, 'kmean = 5 and ', window_size, ' bp windows ', sep = '')

png('hsap_sample_mmus_sample_5cluster.png', width=1200, height=800, units="px")
ggplot(data, aes(x = MDS_1, y = MDS_2, shape = as.factor(factors), colour = as.factor(cluster5$cluster))) + 
  geom_point(size = 2, alpha = 1) +
  labs(title=plot_title, x ="Coordinate 1", y = "Coordinate 2") +
  theme(
    plot.title = element_text(size = 20, face="bold"),
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 15),
    legend.title = element_text(size=18, face = 'bold'),
    legend.text = element_text(size = 15),
    legend.direction = "horizontal",
    legend.position="bottom",
    legend.box = "vertical"
  ) +
  scale_shape_discrete(name = 'Factors', labels = levels(as.factor(factors))) +
  scale_colour_manual(name = 'Clustering', values = rev(colours)) 
dev.off()

png('hsap_sample_mmus_sample_factors.png', width=1200, height=800, units="px")
ggplot(data, aes(x = MDS_1, y = MDS_2, shape = as.factor(species), colour = as.factor(factors))) + 
  geom_point(size = 2, alpha = 1) +
  labs(title=plot_title, x ="Coordinate 1", y = "Coordinate 2") +
  theme(
    plot.title = element_text(size = 20, face="bold"),
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 15),
    legend.title = element_text(size=18, face = 'bold'),
    legend.text = element_text(size = 15),
    legend.direction = "horizontal",
    legend.position="bottom",
    legend.box = "vertical"
  ) +
  scale_colour_manual(name = 'Factors', values = colours, labels = levels(as.factor(factors))) +
  scale_shape_discrete(name = 'Species', labels = species_nice_names) 
dev.off()
