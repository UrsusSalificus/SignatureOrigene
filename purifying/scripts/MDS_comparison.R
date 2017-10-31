#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(ggplot2)

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
all_names <- row.names(fit$points)
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

plot_title <- paste('Distance matrices of either pure factor sequences or whole genome,\nwith kmer = ', 
                    kmer, ' and ', window_size, ' bp windows ', sep = '')

png(output, width=900, height=650, units="px")
ggplot(data, aes(x = MDS_1, y = MDS_2, shape = as.factor(species), colour = as.factor(cut))) + 
  geom_point(size = 2) +
  labs(title=plot_title, x ="Coordinate 1", y = "Coordinate 2") +
  theme(
    plot.title = element_text(size = 20, face="bold"),
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 15),
    legend.title = element_text(size=18, face = 'bold'),
    legend.text = element_text(size = 15)
  )+
  scale_colour_manual(name = 'Factors', values = colours, labels = levels(as.factor(factors))) +
  scale_shape_discrete(name = 'Species', labels = c('plop', 'plap'))
dev.off() 

# TESTING PART
cluster <- hclust(distance_matrix, method = 'ward.D2')
cut <- cutree(cluster, k = 5)

table(species, cut)
table(factors, cut)

ggplot(data, aes(x = MDS_1, y = MDS_2, shape = as.factor(species), colour = as.factor(cut))) + 
  geom_point(size = 2, alpha = 0.1) +
  labs(title=plot_title, x ="Coordinate 1", y = "Coordinate 2") +
  theme(
    plot.title = element_text(size = 20, face="bold"),
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 15),
    legend.title = element_text(size=18, face = 'bold'),
    legend.text = element_text(size = 15)
  )+
  scale_shape_discrete(name = 'Species', labels = c('plop', 'plap'))
