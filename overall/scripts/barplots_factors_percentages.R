#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(ggplot2)
library(ggrepel)

### Comparing the percentages of the different features 
# Will input the arguments:
# 1. Output file
# 2. Window size
args <- commandArgs(trailingOnly=TRUE)

output <- args[1]
window_size <- args[2] 

dict <- readRDS('../input/dictionary.RData')
factor_colours <-  readRDS('../input/factor_colours.RData')

# Get all the perentages file
percentage_directory <- paste('../files/factor_percentages/', as.character(window_size), sep = '')
percentages_files <- Sys.glob(paste(percentage_directory, '*', sep = '/'))

# Merge all the matrices into one
percentages <- do.call(rbind, lapply(percentages_files, function(each_file) {
  species_perc <- read.delim(each_file, header = FALSE)
  colnames(species_perc) <- c('species', 'factors', 'perc', 'l_window', 'n_window')
  
  # Keep only the factor bigger than 0.1%
  species_perc <- species_perc[species_perc$perc > 0.1 & species_perc$n_window > 0,]
  
  # For visual ease, we change very low % in 0 -> we may thus not have perfectly 100
  sum_perc <- sum(species_perc$perc)
  if (sum_perc != 100) {
    uncat <- species_perc$perc[species_perc$factors == 'uncategorized']
    species_perc$perc[species_perc$factors == 'uncategorized'] <- uncat + (100 - sum_perc)
  }
  
  # We will also need to know here to put the labels
  order_by_factor <- species_perc[order(species_perc$factors, decreasing = TRUE),]
  y_for_label <- vector(mode = 'numeric', length = nrow(species_perc))
  each_level <- 1
  to_add <- 0
  while (each_level < nrow(species_perc) + 1) {
    y_for_label[each_level] <- (order_by_factor$perc[each_level] / 2) + to_add
    to_add <- to_add + order_by_factor$perc[each_level]
    each_level <- each_level + 1
  }
  order_by_factor$label <- y_for_label
  
  return(order_by_factor)
}))
colnames(percentages) <- c('species', 'factors', 'perc', 'l_window', 'n_window', 'y_labels')

# We will also prepare the colours
factors_names <- as.character(sort(unique(percentages$factors)))
all_colours = unlist(sapply(factors_names, function (each_factor) {
  # First case, if it is individual factor, find their colours
  if (length(strsplit(each_factor, '-')[[1]]) == 1){
    return(toString(factor_colours$colours[factor_colours$factors == each_factor]))
  }
  # Finally if it is more than one factor, must find compromise between the colours
  else {
    first_factor_col <- toString(factor_colours$colours[factor_colours$factors == strsplit(each_factor, '-')[[1]][1]])
    second_factor_col <- toString(factor_colours$colours[factor_colours$factors == strsplit(each_factor, '-')[[1]][2]])
    # Take the intermediate colour
    col <- colorRampPalette(c(first_factor_col, second_factor_col))(3)[2]
    i <- 2
    # Redo this for each remaining associating factor
    while (i < length(strsplit(each_factor, '-')[[1]])) {
      i_factor_col <- toString(factor_colours$colours[factor_colours$factors == strsplit(each_factor, '-')[[1]][i]])
      col <- colorRampPalette(c(col, i_factor_col))(3)[2]
      i <- i + 1
    }
    return(col)
  }
}))

plot_title <- 'Factors percentages within the whole genome'

# We will order the x axis by genome size
size_genome <- sapply(levels(percentages$species), function(each) {sum(percentages$n_window[percentages$species == each])})
abbrev_right_species <- levels(percentages$species)[order(size_genome)]
percentages$species <- factor(percentages$species, levels = abbrev_right_species)

# Now find the nice names
levels(percentages$species) <- unname(sapply(levels(percentages$species), function (each) {
  return(dict$true[dict$abbrev == each])
}))

# Finally, will also keep only one decimal for visiblity
percentages$perc <- round(percentages$perc, 1) 

png(output, width=1600, height=900, units="px")
ggplot(percentages, aes(x = species, y = perc, fill = factors, label = perc)) +
  geom_bar(stat = "identity") +
  geom_label_repel(aes(label = perc, x = species, y = y_labels), show.legend = FALSE, size = 6) +
  labs(title=plot_title, x ="Species", y = "Percentage [%]") +
  theme(
    plot.title = element_text(size = 20, face="bold"),
    axis.title.x = element_text(size = 20),
    axis.text.x  = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    axis.text.y  = element_text(size = 18),
    legend.title = element_text(size=19, face = 'bold'),
    legend.text = element_text(size=18),
    legend.key = element_blank()
  ) +
  scale_fill_manual(name = 'Factors', values = all_colours, labels = factors_names) +
  guides(fill = guide_legend(ncol = 1))
dev.off() 
