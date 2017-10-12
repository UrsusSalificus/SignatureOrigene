#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(ggplot2)

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

distance_directory <- paste("files/distances/pearson/", window_size, '_', kmer, '/', sep ='')

masked_distance_matrices_files <- Sys.glob(paste(distance_directory, species, '*_masked_*', sep = ''))
pure_distance_matrices_files <- Sys.glob(paste(distance_directory, species, '*_pure_*', sep = ''))
whole_distance_matrix_file <- paste(distance_directory, species, '_whole_dist_matrix.RData', sep = '')

get_mean_distance_masked_or_pure <- function (distance_matrix_path) {
  distance_matrix <- as.matrix(readRDS(distance_matrix_path))
  # We will have to take only the row and column of the factor (and not of the center)
  record_ids <- row.names(distance_matrix)
  i = 1
  # While we have a row name which contain the factor inside it (factor_NC_*somenumbers*), we have row of factor
  while (length(strsplit(record_ids[i], '_')[[1]]) != 2){
    i <- i + 1
  }
  # We have to take the last +1, as it triggered the end of the while
  i <- i -1 
  # We pick only the factor rows, and only the center column ot have the distance of factor VS center only
  # (and not the mean distance among factor or center)
  factor_only_matrix <- distance_matrix[1:i,i:ncol(distance_matrix)]
  # Finally, we can get the raw mean of distance to center of this factor
  # Note: if we have only one "row" (one window), we can directly compute the mean of it
  if (i == 1) {row_mean <- mean(factor_only_matrix)} else {row_mean <- rowMeans(factor_only_matrix)}
  return(row_mean)
} 

get_mean_distance_whole <- function (distance_matrix_path) {
  distance_matrix <- as.matrix(readRDS(Sys.readlink(distance_matrix_path)))
  mean_dist <- rowMeans(distance_matrix)
  sorted <- order(mean_dist)
  center <- sorted[1:round(length(sorted) * 0.1)]
  # This matrix will contain all the row of the windows while having only the column of the center
  # This gives us the distance of all windows (including center) to center
  centered_matrix <- distance_matrix[,center]
  
  rowMeans(centered_matrix)
}  

masked_mean_distance <- sapply(masked_distance_matrices_files, get_mean_distance_masked_or_pure)
pure_mean_distance <- sapply(pure_distance_matrices_files, get_mean_distance_masked_or_pure)
whole_mean_distance <- get_mean_distance_whole(whole_distance_matrix_file)

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
all_mean_factors <- as.factor(c(masked_factors, pure_factors, whole_factors))

# Finally, we need to know if masked or pure
type <- c(rep('Masked', time = length(unlist(masked_mean_distance))),
          rep('Pure', time = length(unlist(pure_mean_distance))),
          rep('NA', time = length(unlist(pure_mean_distance))))



distance_matrix <- as.matrix(distance_matrix)

ncol(distance_matrix)

colnames(distance_matrix)[81]
factor_only_matrix <- distance_matrix[81:ncol(distance_matrix), 1:81]

row_mean <- rowMeans(factor_only_matrix)

all_mean_distance <- c(unlist(masked_mean_distance), unlist(pure_mean_distance), row_mean)
whole_factors <- rep('Whole Genome', time = length(row_mean))
all_mean_factors <- as.factor(c(masked_factors, pure_factors, whole_factors))





dat <- data.frame(mean_dist = all_mean_distance, factors = all_mean_factors)

length(dat$mean_dist[dat$factors == 'CDS'])



# We want the boxplot to be in median order, except for the far left which will always be Whole genome
factor_median <- by(dat$mean_dist, dat$factors, median)[-nlevels(dat$factors)]
factors_ordered <- order(factor_median)
levels_ordered <- names(factor_median)[factors_ordered]
dat$factors <- factor(dat$factors, levels = c('Whole genome', levels_ordered))

# We will also prepare the colours
all_colours <- c("#53b125", "#0057d2", "#d99b24", "#d940ad", "#cf290e")
# We want the whole genome to be black, so we pick the number of non-whole colour, then add black at the end
colours <- c('grey', all_colours[factors_ordered])

plot_title <- paste('Distance to center when masking factor,\nwith kmer = ', 
                    kmer, ' and ', window_size, ' bp windows ', sep = '')





#png(output, width=900, height=650, units="px")
ggplot(dat) + 
  geom_violin(aes(x = factors, y = mean_dist, fill = factors), trim = FALSE) +
  geom_boxplot(aes(x = factors, y = mean_dist), width = 0.1, outlier.colour = 'NA') +
  labs(title=plot_title, x ="Masked factor", y = "Distance to center") +
  theme(
    plot.title = element_text(size = 20, face="bold"),
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 15)
  )
#+ geom_jitter(aes(x = factors, y = mean_dist), shape=16, position=position_jitter(0.3), alpha = 0.1)
#dev.off() 

