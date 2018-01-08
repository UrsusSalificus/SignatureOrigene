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
# 6. sample size
# 7. species name
# 8. comparison species name
args <- commandArgs(trailingOnly=TRUE)

output = args[1]
window_size <- args[4]
kmer <- args[5]
sample_size <- args[6]
species_abbrev <- args[7]
comparison_abbrev <- args[8]

dict <- readRDS('../input/dictionary.RData')
factor_colours <-  readRDS('../input/factor_colours.RData')

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
factors_names <- sort(unique(factors))


# We will also prepare the colours
all_colours = unlist(sapply(factors_names, function (each_factor) {
  # First case, if factor is Whole Genome -> black
  if (each_factor == 'whole') {return('black')}
  # Else, if it is individual factor, find their colours
  else if (length(strsplit(each_factor, '-')[[1]]) == 1){
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

data <- data.frame(MDS_1 = fit$points[,1], MDS_2 = fit$points[,2], factors = factors, species = species)

species_nice_names <- unname(sapply(levels(as.factor(species)), function (each) {
  return(dict$true[dict$abbrev == each])
}))

plot_title <- paste('Distance between ', species_nice_names[1], '/', species_nice_names[2], 
                    ' pure feature sequences,\nWith maximum ', sample_size, ' windows for each feature, kmer = ', 
                    kmer, ' and ', window_size, ' bp windows', sep = '')

png(output, width=30, height=18, units="cm", res = 300)
ggplot(data, aes(x = MDS_1, y = MDS_2, shape = as.factor(species), fill = as.factor(factors))) + 
  geom_point(size = 2, alpha = 1) +
  labs(x ="Coordinate 1", y = "Coordinate 2") +
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 15),
    legend.title = element_text(size=18, face = 'bold'),
    legend.text = element_text(size = 15)
  ) +
  guides(
    fill = guide_legend(ncol = 1, override.aes = list(shape = 21)),
    shape = guide_legend(override.aes = list(fill = 'black'))
  ) +
  scale_shape_manual(name = 'Species', values = c(21, 24), labels = species_nice_names) +
  scale_fill_manual(name = 'Features', values = all_colours, labels = factors_names)
dev.off() 


