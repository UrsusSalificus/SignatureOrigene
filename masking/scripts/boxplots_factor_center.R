#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(ggplot2)
library(RColorBrewer)

### Boxplot comparison of mean distance to center
# Will input the arguments:
# 1. path to the output image (with its name in the path)
# 2. window size
# 3. k-mer size
# 4. species 
args <- commandArgs(trailingOnly=TRUE)

output = args[1]
window_size <- args[2]
kmer <- args[3]
species <- args[4]

# We will need a dictionnary of species abbreviation VS nice species name
dict <- readRDS('../input/dictionary.RData')

distance_directory <- paste("files/distances/pearson/", window_size, '_', kmer, '/', sep ='')

# Matrix files' path
masked_distance_matrices_files <- Sys.glob(paste(distance_directory, species, '*_masked_*', sep = ''))
pure_distance_matrices_files <- Sys.glob(paste(distance_directory, species, '*_pure_*', sep = ''))
whole_distance_matrix_file <- paste(distance_directory, species, '_whole_vs_center_dist_matrix.RData', sep = '')

get_mean_distance_masked_or_pure <- function (distance_matrix_path) {
  distance_matrix <- as.matrix(readRDS(distance_matrix_path))
  # We will have to take only the row and column of the factor (and not of the center)
  # To do this we will exculde the center rows/columns 
  center <- which(row.names(distance_matrix) == 'center')
  
  # This matrix will contain all the row of the factor while having only the column of the center
  factor_only_matrix <- distance_matrix[-center, center]
  # Finally, we can get the raw mean of distance to center of this factor
  # Note: if we have only one "row" (one window), we can directly compute the mean of it
  if (ncol(distance_matrix) - length(center) == 1) {
    row_mean <- median(factor_only_matrix)
  } else {
    row_mean <- apply(factor_only_matrix, 1, median)
  }
  return(row_mean)
} 

masked_mean_distance <- sapply(masked_distance_matrices_files, get_mean_distance_masked_or_pure)
pure_mean_distance <- sapply(pure_distance_matrices_files, get_mean_distance_masked_or_pure)
whole_mean_distance <- get_mean_distance_masked_or_pure(whole_distance_matrix_file)

# We will concatenate all the distance for the data frame afterward
all_mean_distance <- c(unlist(masked_mean_distance), unlist(pure_mean_distance), whole_mean_distance)

# We will now create a factor vector:
get_factor_masked_or_pure <- function(distance_matrix_files, mean_distance, matrix_number) {
  file <- distance_matrix_files[matrix_number]
  factor <- strsplit(file, '_')[[1]][4]
  factor_as_vector <- rep (factor, time = length(mean_distance[[matrix_number]]))
  return(factor_as_vector)
}

masked_factors <- unlist(sapply(1:length(masked_distance_matrices_files), function(each_matrix) {
  get_factor_masked_or_pure (masked_distance_matrices_files, masked_mean_distance, each_matrix)
}))
pure_factors <- unlist(sapply(1:length(pure_distance_matrices_files), function(each_matrix) {
  get_factor_masked_or_pure (pure_distance_matrices_files, pure_mean_distance, each_matrix)
}))
whole_factors <- rep('Whole Genome', time = length(whole_mean_distance))

# Same for the factors:
all_factors <- as.factor(c(masked_factors, pure_factors, whole_factors))

# We also need to know if masked or pure
type <- c(rep('Masked', time = length(unlist(masked_mean_distance))),
          rep('Pure', time = length(unlist(pure_mean_distance))),
          rep('NA', time = length(whole_factors)))

# Finally, we need a merged factor of the two:
merged <- factor(sapply(1:length(all_factors), function(each_factor){
  paste(all_factors[each_factor], type[each_factor])
}))

dat <- data.frame(mean_dist = all_mean_distance, factors = all_factors, type = type, merged = merged)

plot_title <- paste(dict$true[dict$abbrev == species], ': distance to center when masking factor\nWith kmer = ', 
                    kmer, ' and ', window_size, ' bp windows ', sep = '')

# If we want to label the number of windows per group:
# Write function which will be computed on each group using ggplot2 stat_summary
give.n <- function(x){
  return(c(y = max(x) * 1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

png(output, width=900, height=650, units="px")
ggplot(dat,aes(x = factors, y = mean_dist, grp = merged)) + 
  # Invisible points -> enable us to control the color legend to have legend we want
  # Note the shape 15 -> square for square in legend...
  geom_point(aes(col = type), alpha = 0, shape = 15, size = 6) +
  geom_violin(aes(fill = type), trim = FALSE) +
  geom_boxplot(aes(fill = type), width = 0.2, outlier.colour = 'NA', position=position_dodge(0.9)) +
  labs(title=plot_title, x ="Factors", y = "Distance to center") +
  theme(
    plot.title = element_text(size = 18, face="bold"),
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 15),
    legend.title = element_text(size=16, face = 'bold'),
    legend.text = element_text(size=15),
    legend.key = element_blank()    # Will remove the grey background
  ) +
  # Change the visible filling colour
  scale_fill_manual(values = c('#EE6B6B','grey', '#537DE6')) + 
  # Change the invisible (the legend boxes) colours and labels (note that can avoid the Whole grey color label)
  scale_color_manual(values = c('#EE6B6B','grey', '#537DE6'), breaks=c("Masked","Pure")) + 
  guides(grp = FALSE) +
  # Remove the fill legend
  guides(fill = FALSE) + 
  # Legend we want
  guides(col = guide_legend(override.aes = list(alpha = 1), title = 'Type of\nsequence')) +
  stat_summary(fun.data = give.n, geom = "text", position = position_dodge(width = 0.9))
dev.off() 





      

