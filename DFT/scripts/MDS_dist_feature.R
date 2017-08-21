#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(RColorBrewer)

### MultiDimensional Scaling (MDS) analysis of a distance matrix
# Will input the arguments:
# 1. path to the distance matrix
# 2. path to the feature file
# 3. feature as displayed in the title and y axis
# 4. distance matrix information (type of genomic signature/type of distance/window size)
# 5. path to the output image (with its name in the path)
args <- commandArgs(trailingOnly=TRUE)

feature_type = basename(args[3])
window_size = basename(args[4])
output = args[6]

distance_matrix <- read.table(args[1], sep = "\t", header = FALSE)
# Will remove any empty column
distance_matrix <- distance_matrix[, colSums(is.na(distance_matrix)) == 0]

# Copying the upper diagonal into the lower one
distance_matrix[lower.tri(distance_matrix)] <- t(distance_matrix)[lower.tri(distance_matrix)]

# Finding which column (= which region) as the smallest distance to the other column/region.
# This column/region will then represent the median region.
min_column = which.min(colSums(distance_matrix))
median_region = distance_matrix[-min_column, min_column]

# MDS
fit <- cmdscale(distance_matrix,eig=TRUE, k=2) # k is the number of dim


feature = read.csv(args[2], sep = '\t', header = FALSE)

# We will check how much we can round the result and still not lose significant differences:
size <- feature[,2]
good_decimal <- 1
round_up <- TRUE
while (round(min(size), good_decimal) == 0) {
  good_decimal = good_decimal + 1
  # If the minimum region contain less than 0.1% CDS, stop here
  if (good_decimal == 3){
    break
  }
}
size <- round(size, good_decimal)

# Colours: if bigger than 10 -> must chose special brewer palette
if (length(size) >= 10) {
  # Will add more feature than needed, to make sure
  colours_heat = rep(rev(brewer.pal(10, 'RdYlBu')), each=ceiling(length(size)/10))
} else {
  colours_heat = rev(brewer.pal(size, 'RdYlBu'))
}
increasing_order = order(size)

plot_title <- paste('MDS metrics of DFTs/', window_size, ' bp windows distance matrix', sep = '')

# Will just plot the result to get good x limits
plot(fit$points[,1], fit$points[,2])
range_x <- par("usr")[1:2][2] - par("usr")[1:2][1]
if (min(fit$points[,1]) < 0) {
  good_x_lim <- c( min(fit$points[,1]) - range_x * 0.25, par("usr")[1:2][2])
} else {
  good_x_lim <- c( min(fit$points[,1]) + range_x * 0.25, par("usr")[1:2][2])
}

# Looking for a good scales of number for the legend:
scale <- vector("integer", length = 10)
for (each_step in 1:10) {
  scale[each_step] <- paste(size[increasing_order][floor(length(size)/10) * each_step], "% of", feature_type)
}
scale <- rev(scale)

png(output, width=700, height=500, units="px")
plot(fit$points[,1], fit$points[,2], xlab="Coordinate 1", ylab="Coordinate 2", main=plot_title, xlim=good_x_lim,
     pch=21, col = 'black', bg = colours_heat[increasing_order], cex =1.5)
legend('topleft',  scale, pch =21, pt.bg = brewer.pal(10, 'RdYlBu'), pt.cex=1.5)
dev.off() 

# Not sure why, it also save the good x limit plot as pdf...
file.remove("Rplots.pdf")


distance_matrix <- matrix(0, ncol = 10, nrow = 10)
distance_matrix[upper.tri(distance_matrix)] <- c(1,2,4,5,6,7, sample(95:100, 39, replace = TRUE))
feature <- data.frame(rep(0, time = 10), c(sample(90:100, 4, replace = TRUE), sample(2:10, 6, replace = TRUE)))


