#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

suppressMessages(library(ggplot2))
suppressMessages(library(dunn.test))
suppressMessages(library(igraph))
library(grid)

### Boxplot comparison of mean distance to center
# Will input the arguments:
# 1. path to the output image (with its name in the path)
# 2. window size
# 3. k-mer size
# 4. species 
# 5. number of samples
args <- commandArgs(trailingOnly=TRUE)

output = args[1]
window_size <- args[2]
kmer <- args[3]
species <- args[4]
n_samples <- args[5]

# We will need a dictionnary of species abbreviation VS nice species name
dict <- readRDS('../input/dictionary.RData')

distance_directory <- paste("files/distances/manhattan/", window_size, '_', n_samples, '_', kmer, sep ='')

# Method of p.values adjustment
method_pval = 'bh'



#### Processing distance matrices ####

# Matrix files' path
masked_distance_matrices_files <- Sys.glob(paste(distance_directory, species, '*_masked_*', sep = '/'))
pure_distance_matrices_files <- Sys.glob(paste(distance_directory, species, '*_pure_*', sep = '/'))
whole_distance_matrix_file <- paste(distance_directory, species, 'whole_vs_center_dist_matrix.RData', sep = '/')

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
    names(row_mean) <- colnames(distance_matrix)[1]
  } else {
    row_mean <- apply(factor_only_matrix, 1, median)
  }
  return(row_mean)
} 

masked_mean_distance <- lapply(masked_distance_matrices_files, get_mean_distance_masked_or_pure)
pure_mean_distance <- lapply(pure_distance_matrices_files, get_mean_distance_masked_or_pure)
whole_mean_distance <- get_mean_distance_masked_or_pure(whole_distance_matrix_file)

# We will concatenate all the distance for the data frame afterward
all_mean_distance <- c(unlist(masked_mean_distance), unlist(pure_mean_distance), whole_mean_distance)

# We will now create a factor vector:
masked_factors <- unlist(lapply(masked_mean_distance, names))
pure_factors <- unlist(lapply(pure_mean_distance, names))
whole_factors <- rep('Whole Genome', time = length(whole_mean_distance))

# Make them as factor
all_factors <- as.factor(c(masked_factors, pure_factors, whole_factors))

# We also need to know if masked or pure
type <- c(rep('Absent', time = length(unlist(masked_mean_distance))),
          rep('Only', time = length(unlist(pure_mean_distance))),
          rep('NA', time = length(whole_factors)))

# Finally, we need a merged factor of the two:
merged <- factor(sapply(1:length(all_factors), function(each_factor){
  paste(all_factors[each_factor], type[each_factor])
}))

dat <- data.frame(mean_dist = all_mean_distance, factors = all_factors, type = type, merged = merged)

plot_title <- paste(dict$true[dict$abbrev == species], 
                    ": impact of factor's presence/absence on the distance to center\nWith kmer = ", 
                    kmer, ' and maximum ', n_samples, ' -',window_size, ' bp long- samples', sep = '')

# If we want to label the number of windows per group:
# Write function which will be computed on each group using ggplot2 stat_summary
give.n <- function(x){
  return(c(y = max(x) * 1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}



#### Statistical tests ####

# This function will extract the groups which are similar following a post-hoc statistical test
grouping <- function (p_values, factors, response) {
  #Function built with help found on 
  #http://menugget.blogspot.ch/2014/05/automated-determination-of-distribution.html#more
  # First find the nuzmber of factors to build the comparison matrix
  n_tot <- nlevels(factors)
  
  # This matrix will compare all factor -> if they are significantly different -> 0, else 1
  similar_matrix <- matrix(nrow=n_tot,ncol=n_tot)
  similar_matrix[upper.tri(similar_matrix)] <- p_values > 0.05
  similar_matrix <- replace(similar_matrix, is.na(similar_matrix), FALSE) # replace NAs with FALSE
  similar_matrix <- similar_matrix + t(similar_matrix)
  diag(similar_matrix) <- 1
  
  # As there will be different versions of this matrix, rename it g for ease
  g <- similar_matrix
  rownames(g) <- 1:n_tot # change row names to assigned number
  colnames(g) <- 1:n_tot # change column names to assigned number
  # Find all features which are similar
  same <- which(g==1)
  # Build a matrix which will store all the different similar couples (e.g. 1-2, 1-3, 3-4, ...)
  g2 <- data.frame(N1=((same-1) %% n_tot) + 1, N2=((same-1) %/% n_tot) + 1)
  g2 <- g2[order(g2[[1]]),] # Get rid of loops (e.g. 1-2 and 2-1) and ensure right naming of vertices
  # Build a graph of all the similar vertices 
  g3 <- simplify(graph.data.frame(g2,directed = FALSE))
  
  # Find all the cliques -> all the vertices which form a distinct group from others
  cliq <- maximal.cliques(g3)
  # We also sort it so the smaller the letter the smaller the distance
  cliq_mean <- lapply(cliq, as.numeric)
  med <- as.numeric(unlist(by(response, factors, function(x)  median (x,na.rm=TRUE))))
  calc_grp_med <- function (x) mean(med[x])
  order_grp <- as.numeric(lapply(cliq_mean, function(x) calc_grp_med(x)))
  cliq <- cliq_mean[order(order_grp,decreasing=FALSE)]
  
  groups <- vector(mode="list", n_tot) # empty list
  lab <- letters[seq(cliq)] # clique labels
  for(i in seq(cliq)){ # loop to concatenate clique labels
    for(j in cliq[[i]]){
      groups[[j]] <- paste0(groups[[j]], lab[i])
    }
  }
  groups <- as.factor(unlist(groups))
  
  return(groups)
}

comp_features <- data.frame(features = droplevels(as.factor(dat$factors[dat$type == 'Only'])), 
                            distance = as.numeric(dat$mean_dist[dat$type == 'Only']))

# We will write the dunn test in between features on a file
comp_feature_file <- paste(dirname(output), '/', species, '_dunn_test_between_features.txt', sep = '')
# Note, we write the default test, without p_value adjustment
write(capture.output(dunn.test(x = comp_features$distance, g = comp_features$features)), comp_feature_file)

# Benjamin-Hochberg adjusted p-values of a two-tailed (*2) Dunn.test
p_values <- dunn.test(x = comp_features$distance, g = comp_features$features, list = FALSE,
                      table = FALSE, kw = FALSE, method = method_pval)$P.adjusted * 2
groups <- grouping(p_values = p_values, factors = comp_features$features, response = comp_features$distance)
# Add the empty Whole genome
groups <- as.character(groups) ; groups[length(groups) + 1] <- ''

# Then we will also perform absent/only comparisons (note we remove the whole genome)
comp_state <- split(dat, dat$factors) ; comp_state <- comp_state[1:length(comp_state) - 1]
# Same for wilcoxon test -> store them in a file
comp_state_file <- paste(dirname(output), '/', species, '_wilcoxon_each_feature.txt', sep = '')
write(capture.output(lapply(comp_state, function(each) wilcox.test(each$mean_dist~each$type))), file = comp_state_file)
# We extract the p.values to add significant stars to the plot
state_p_values <- unlist(lapply(comp_state, function(each) wilcox.test(each$mean_dist~each$type)$p.val))

# Assigning stars to p.values 
star_check <- function(p_val) {
  if (p_val < 0.05 & p_val > 0.001) {
    return('*')
  } else if (p_val < 0.001 & p_val > 0.0001) {
    return('**')
  } else if (p_val < 0.0001) {
    return('***')
  } else {
    return('')
  }
}

stars_comp_state <- sapply(state_p_values, star_check)
# Add the empty Whole genome
stars_comp_state[length(stars_comp_state) + 1] <- ''

# This plot will help us know how much we must add to get to the top of the boxplot
p <- ggplot(dat,aes(x = factors, y = mean_dist, grp = merged)) +
  geom_violin(aes(fill = type), trim = FALSE)

y_ranges <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range
in_between <- abs(y_ranges[2]-y_ranges[1])
to_add <- in_between*0.1

# And we also need the max/min distance for each feature
max_each_feature = sapply(split(dat$mean_dist, dat$factors), max)



#### Plotting ####

png(output, width=900, height=650, units="px")
ggplot(dat,aes(x = factors, y = mean_dist, grp = merged)) + 
  # Invisible points -> enable us to control the color legend to have legend we want
  # Note the shape 15 -> square for square in legend...
  geom_point(aes(col = type), alpha = 0, shape = 15, size = 5) +
  geom_violin(aes(fill = type), trim = FALSE) +
  # First annotation, at the bottom -> the stars
  annotate('text', x = 1:nlevels(dat$factors), y = max_each_feature + to_add, angle = 0, 
           label = stars_comp_state, size = 6,  fontface="bold") +
  annotate('text', x = 1:nlevels(dat$factors), y = max_each_feature + to_add * 1.5, angle = 0, 
           label = groups, size = 6,  fontface="bold") +
  geom_boxplot(aes(fill = type), width = 0.2, outlier.colour = 'NA', position=position_dodge(0.9)) +
  labs(title=plot_title, x ="Factors", y = "Distance to center") +
  theme(
    plot.title = element_text(size = 18, face="bold"),
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 15, angle = 90, hjust = 1),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 15),
    legend.title = element_text(size=16, face = 'bold'),
    legend.text = element_text(size=15),
    legend.key = element_blank()    # Will remove the grey background
  ) +
  # Change the visible filling colour
  scale_fill_manual(values = c('#537DE6','grey', '#EE6B6B')) + 
  # Change the invisible (the legend boxes) colours and labels (note that can avoid the Whole grey color label)
  scale_color_manual(values = c('#537DE6','grey', '#EE6B6B'), breaks=c("Absent","Only")) + 
  guides(grp = FALSE) +
  # Remove the fill legend
  guides(fill = FALSE) + 
  # Legend we want
  guides(col = guide_legend(override.aes = list(alpha = 1), title = "Factor's state\nin the sequences")) +
  stat_summary(fun.data = give.n, geom = "text", position = position_dodge(width = 0.9))
dev.off() 
