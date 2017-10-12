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

distance_matrices <- Sys.glob(paste(distance_directory, species, '*', sep = ''))
mean_distances <- lapply(distance_matrices, function(each_matrix) {
  distance_matrix <- as.matrix(readRDS(each_matrix))
  if (strsplit(each_matrix, '_')[[1]][4] != 'whole'){
    # We will have to take only the row and column of the factor (and not of the center)
    record_ids <- row.names(distance_matrix)
    i = 1
    # While we have a row name which contain the factor inside it (factor_NC_*somenumbers*), we have row of factor
    while (length(strsplit(record_ids[i], '_')[[1]]) != 2){
      i <- i + 1
    }
    # We have to take the last +1, as it triggered the end of the while
    i <- i -1 
    # As it is in the same order, we can use it for column as well
    factor_only_matrix <- distance_matrix[1:i,1:i]
    # Finally, we can get the raw mean of distance to center of this factor
    rowMeans(factor_only_matrix)
  } else {
    # Else, we are working on the whole genome file, and we must find the center, and find how much the other 
    # windows are different from it
    mean_dist <- rowMeans(distance_matrix)
    sorted <- order(mean_dist)
    center <- sorted[1:round(length(sorted) * 0.1)]
    # This matrix will contain all the row of the windows outside of the center
    # while having only the column of the center
    # This gives us the distance of all windows to center
    centered_matrix <- distance_matrix[-center,center]
    
    rowMeans(centered_matrix)
  }
})

# We will extract the factor anme through the file name and length of row mean of each matrix
factors <- unlist(sapply(1:length(distance_matrices), function(each_matrix) {
  file <- distance_matrices[each_matrix]
  if (strsplit(file, '_')[[1]][4] != 'whole'){
    factor <- strsplit(file, '_')[[1]][4]
    rep (factor, time = length(mean_distances[[each_matrix]]))
  } else {
    rep('Whole genome', time = length(mean_distances[[each_matrix]]))
  }
}))

dat <- data.frame(mean_dist = unlist(mean_distances), factors = as.factor(factors))

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

png(output, width=900, height=650, units="px")
ggplot(dat) + 
  geom_violin(aes(x = factors, y = mean_dist, fill = factors), trim = FALSE) +
  scale_fill_manual(values = colours, guide=FALSE) +
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
dev.off() 

