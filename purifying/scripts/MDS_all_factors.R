#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(ggplot2)

### MultiDimensional Scaling (MDS) analysis of a distance matrix
# Will input the arguments:
# 1. path to the output image (with its name in the path)
# 2. path to the fit of the distance matrix
# 3. window size
# 4. k-mer size
# 5. species name
# 6. sample size
args <- commandArgs(trailingOnly=TRUE)

output = args[1]
window_size <- args[3]
kmer <- args[4]
species <- args[5]
sample_size <- args[6]

# We will need a dictionnary of species abbreviation VS nice species name
dict <- readRDS('../input/dictionary.RData')
factor_colours <-  readRDS('../input/factor_colours.RData')

# MDS
fit <- readRDS(args[2])

# We now have to fetch from which factor it is from
all_names <- row.names(fit$points)
factors <- unname(sapply(all_names, function(each_line) {
  if (length(strsplit(each_line, '_')[[1]]) == 3) {
    # If when split we get 3 result = untouched genome (only the record name) + start of the window
    return('whole')
  } else {
    # Else, it's the pure samples, hence colorRampPalettefactor abbreviation -> keep those
    return(each_line)
  }
}))
factors_names <- unique(factors)

# We will also prepare the colours
# The first vector will contain the "basic" colors -> one for each factor
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

# We will reduce visibility of the whole genome points:
reduced_vis <- unname(sapply(factors, function (each_line) {
  if (each_line == 'whole'){
    return(TRUE)
  } else {
    return(FALSE)
  }
}))

data <- data.frame(MDS_1 = fit$points[,1], MDS_2 = fit$points[,2], factors = factors, reduced_vis = reduced_vis)

plot_title <- paste(dict$true[dict$abbrev == species], 
                    ': proximity between the whole genome and feature only sequences\nWith maximum ', 
                    sample_size, ' windows for each feature, kmer = ', kmer, ' and ', window_size, ' bp windows ', sep = '')

png(output, width=30, height=18, units="cm", res = 300)
ggplot(data, aes(x = MDS_1, y = MDS_2, colour = as.factor(factors), alpha = reduced_vis)) + 
  geom_point(size = 2) +
  labs(x ="Coordinate 1", y = "Coordinate 2") +
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 15),
    legend.title = element_text(size=18, face = 'bold'),
    legend.text = element_text(size = 15)
  ) +
  scale_colour_manual(name = 'Features', values = all_colours,
                      labels = c(levels(as.factor(factors))[-nlevels(as.factor(factors))], 'Whole genome')) +
  scale_alpha_discrete(range = c(1, 0.4), guide=FALSE) + 
  guides(
    colour = guide_legend(ncol = 1)
  )
dev.off() 




