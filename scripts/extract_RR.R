#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(data.table)

### Infer recombination rates, using as input the coordinates files produced by the homonyme Python script.
# Will input the arguments:
# 1. path to the input file = coordinates of each of the windows of the species' genome
# 2. path to the RR spline functions of the species
# 3. path to the output table of recombination rate per windows
args <- commandArgs(trailingOnly=TRUE)

# Load the windows' coordinates
coordinates <- fread(args[1])
# The chromosome is a factor as default, which creates trouble afterward
coordinates$V1 <- as.character(coordinates$V1)

# Load the RR spline fits
load(args[2]) ; fit_spline <- fit_spline

split_coord <- split(coordinates, coordinates$V1)

# We will create a vector of recombination rate as long as there is rows/windows on the coordinates table
windows_rr <- sapply(1:length(split_coord), function(each_chrom) {
  chrom_windows <- split_coord[[each_chrom]]
  chrom_fit <- fit_spline[[which(names(fit_spline) == chrom_windows$V1[1])]]
  # The first and last markers position (bp) usually marks the boundaries of telomeres
  first_marker <- chrom_fit$x[1] ; last_marker <- chrom_fit$x[length(chrom_fit$x)]
  # We want a table at the end, here with sapply, we end up with column = each window, whereas we want the reverse (hence t())
  sapply(1:nrow(chrom_windows), function (each_window) {
    window_line <- chrom_windows[each_window,]
    # If not telomeres, compute average recombination rate
    if (window_line$V3 < last_marker && window_line$V3  > first_marker){
      middle_window <- mean(c(as.numeric(window_line$V2), as.numeric(window_line$V3)))
      predicted_rr <- predict(chrom_fit, middle_window)$y
      # Second level of control: if the recombination rate is < 0, mark as 0 instead
      if (predicted_rr < 0) {
        return(0)
      # Else, can return the recombination rate
      } else {
        return(predicted_rr)
      }
    # Else, in the telomere = 0 recombination rate
    } else {
      return(0)
    }
  })
})

all_windows_rr <- matrix(c(coordinates$V1, unlist(windows_rr)), ncol =2)

write.table(all_windows_rr, args[3], quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
