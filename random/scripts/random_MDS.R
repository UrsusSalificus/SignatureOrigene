#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(data.table)
library(RColorBrewer)
library(ggplot2)

# Importing both species name and abbreviation:
full <- read.table("../input/all_species.txt", sep = "\t", header = FALSE)[-1,1]
abbrev <- unname(unlist(read.table("../input/abbrev_species.txt", sep = " ", header = FALSE)[1,-1]))
translation <- data.frame(full = full, abbrev = abbrev)

# Importing the tables containing all n windows per species FCGRs/DFTs
FCGRs <- fread("temp/FCGRs", sep = "\t", header = FALSE)
DFTs <- fread("temp/DFTs", sep = "\t", header = FALSE)
GS <- list(FCGRs = FCGRs, DFTs = DFTs)

window_size <- ncol(DFTs)-2
kmer <- log(ncol(FCGRs)-2, base = 4)

for (each_gs in 1:length(GS)) {
  gs_table <- as.data.frame(GS[each_gs])
  species <- factor(gs_table[,1])
  n_species <- nlevels(species)
  sample_size <- nrow(gs_table)/n_species
  
  # Transform the abbreviations in full names
  full_species <- ''
  for (each_species in 1:n_species) {
    which_row <- which(as.character(translation$abbrev) == as.character(levels(species)[each_species]))
    full_species[each_species] <- as.character(translation$full[which_row])
  }
  levels(species) <- full_species
 
  # Plot title, output file and distance used varies between type of genomic signature used
  gs_type <- names(GS)[each_gs]
  if (gs_type == 'FCGRs') {
    distance_matrix <- dist(gs_table[,2:(ncol(gs_table)-1)], method = 'manhattan')
    plot_title <- paste("Comparison of", sample_size, "x", window_size, 
                        "bp windows per species\nGenomic signature used:", 
                        gs_type,  "with k-mer = ", kmer)
    output_directory <- paste("files/results/FCGRs/", paste(sample_size, window_size, kmer, sep = '_'), sep = '')
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
    output <- paste(output_directory, "/random_FCGR_MDS.png", sep = '')
  } else {
    distance_matrix <- dist(gs_table[,2:(ncol(gs_table)-1)])
    plot_title <- paste("Comparison of", sample_size, "x", window_size, 
                        "bp windows per species\nGenomic signature used:", gs_type)
    output_directory <- paste("files/results/DFTs/", paste(sample_size, window_size, sep = '_'), sep = '')
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
    output <- paste(output_directory, "/random_DFT_MDS.png", sep = '')
  }
  
  # MDS
  fit <- cmdscale(distance_matrix,eig=TRUE, k=2) # k is the number of dim
  
  # Type of cluster
  cluster <- hclust(distance_matrix, method = "ward.D2")
  cluster_cut <- cutree(cluster, n_species)
  
  #cluster <- kmeans(distance_matrix, n_species)
  #cluster_cut <- cluster$cluster
  
  # Now check if in the right cluster or not
  data <- data.frame(MDS_1 = fit$points[,1], MDS_2 = fit$points[,2], species = species, cluster=cluster_cut)
  # First look for what is the most representative cluster for each species
  species_cluster <- ''
  for (each_species in 1:nlevels(data$species)) {
    count <- table(data$cluster[data$species == unique(data$species)[each_species]])
    species_cluster[each_species] <- names(sort(count, decreasing = TRUE))[1]
  }
  data$true <- rep(as.numeric(species_cluster), each = sample_size)
  right <- apply(data[4:5], 1, function(each_row) {identical(as.integer(each_row[1]), as.integer(each_row[2]))})
  data$right <- unlist(right)
  
  colours <- c("#6839a6","#3a6f2f","#c94e8f","#736221","#006591","#a44023","#8f8ee1")

  # The plot in itself
  png(output, width=700, height=500, units="px")
  plot <- ggplot(data, aes(x = MDS_1, y = MDS_2)) + 
    geom_point(aes(shape = factor(cluster), colour = factor(species)),size = 4) + 
    labs(title=plot_title, x ="Coordinate 1", y = "Coordinate 2") +
    theme(
      plot.title = element_text(size = 15, face="bold"),
      axis.title.x = element_text(size = 15),
      axis.text.x  = element_text(size = 12),
      axis.title.y = element_text(size = 15),
      axis.text.y  = element_text(size = 12)
    ) + scale_shape_manual(name = "Clustering", values=c(0, 1, 2, 3, 7, 8, 9)[1:nlevels(species)]) + 
    scale_color_manual(name = "Species", values=colours[1:nlevels(species)])
  print(plot)
  dev.off() 
}





