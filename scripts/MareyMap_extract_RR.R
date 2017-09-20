#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(MareyMap)

### Extract Recombination rate for A. thaliana, C. elegans and S. cerevisiae
# A. thaliana and C. elegans are already implimanted in the MareyMap package
# S. cerevisiae must be incorporated:
pgmap <- read_delim("tools/RRC/s_cerevisiae/pgmap.xls", "\t")

split_map <- split(pgmap, pgmap$chr)

Saccharomyces_cerevisiae <- new("MapSet")
setName(Saccharomyces_cerevisiae) <- 'Saccharomyces_cerevisiae'

number_chromosome <- 16
# We will not use the 17th (?) chromosome, hence the 1:16
for (each_chrom in 1:16) {
  each_map <- split_map[[each_chrom]]
  # Certain markers do not have genetic distances, we thus get rid of them
  to_remove <- which(is.na(each_map$genetic))
  each_map <- each_map[-to_remove,]
  physical_positions <- apply(each_map, 1, function(each_row) {
    round(mean(as.numeric(each_row[4:5])))
  }) 
  
  # Constructing the MareyMap object
  map_to_add <- MareyMap()
  setName(map_to_add) <- 'Saccharomyces_cerevisiae'
  mapName(map_to_add) <- as.character(each_map$chr[1])
  markerNames(map_to_add) <- as.character(each_map$gene)
  physicalPositions(map_to_add) <- physical_positions
  geneticDistances(map_to_add) <- each_map$genetic
  markerValidity(map_to_add) <- rep(TRUE, nrow(each_map))
  
  # Adding the map to the set
  Saccharomyces_cerevisiae <- Saccharomyces_cerevisiae + map_to_add
}

# 