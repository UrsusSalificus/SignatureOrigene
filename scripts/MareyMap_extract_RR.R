#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(MareyMap)

### Extract spline fitting function of Recombination Rate (RR)

get_chrom_RR <- function(chrom_genetic_map) {
  chrom_genetic_map <- chrom_genetic_map + new("MMSpline3")
  recombination_rates <- chrom_genetic_map@interpolations$spline@rates
  return(recombination_rates)
}

get_all_RR <- function(all_genetic_maps, record_names) {
  # Retrieve the Chromosome IDs (record_names$V1) of the specific species (setName(all_genetic_maps))
  chromosome_ids <- record_names$V1[record_names$V2 == setName(all_genetic_maps)]
  
  # Will append all the different smooth spline fit
  fitted_models <- list()
  
  # Chromosome by chromosome, use MareyMap to find the recombination rates, then fit a smooth spline  
  for (each_chrom in 1:length(all_genetic_maps)) {
    recombination_rates <- get_chrom_RR(all_genetic_maps[[each_chrom]])
    # Any negative rates is most probably due to fitting issues, we will thus change it to no recombination rate
    if (any(recombination_rates < 0)) {
      recombination_rates[recombination_rates < 0] <- 0
    }
    spline <- smooth.spline(all_genetic_maps[[each_chrom]]@physicalPositions, 
                            recombination_rates, cv= TRUE)
    fitted_models[[each_chrom]] <- spline
  }
  # Associate the id to each spline fit
  names(fitted_models) <- chromosome_ids
  return(fitted_models)
}

record_names <- read.table("wanted_records.txt")

# We will also create the directory which will contain all these spline fits
dir.create('data/features')




# - H. sapiens recombination rate was recovered from Hapmap
# - A. thaliana and C. elegans genetic distances data were found in the MareyMap package, 
# and will then be used to compute RRs.

# - S. cerevisiae must be incorporated:



# 1) H. sapiens







# 2) M. musculus
map <- read.csv('tools/RRC/m_musculus/Revised_HSmap_SNPs.csv')

split_map <- split(map, map$chr)

Mus_musculus <- new("MapSet")
setName(Mus_musculus) <- 'Mus_musculus'

number_chromosome <- 20
for (each_chrom in 1:number_chromosome) {
  each_map <- split_map[[each_chrom]]
  
  # Constructing the MareyMap object
  map_to_add <- new("MareyMap")
  setName(map_to_add) <- 'Mus_musculus'
  markerNames(map_to_add) <- as.character(each_map$snpID)
  physicalPositions(map_to_add) <- as.numeric(each_map$build37)
  if (each_chrom != 20) {
    geneticDistances(map_to_add) <- as.numeric(each_map$ave_cM)
    mapName(map_to_add) <- as.character(each_map$chr[1])
  } else {
    geneticDistances(map_to_add) <- as.numeric(each_map$fem_cM)
    mapName(map_to_add) <- 'X'
  }
  markerValidity(map_to_add) <- rep(TRUE, nrow(each_map))
  
  # Adding the map to the set
  Mus_musculus <- Mus_musculus + map_to_add
}

setName(Mus_musculus) <- 'm_musculus'

fit_spline <- get_all_RR(all_genetic_maps = Mus_musculus, record_names = record_names)

save(fit_spline, file = 'data/features/m_musculus_RR_spline.RData')



# 3) C. elegans (already incorporated in MareyMap)
data("Caenorhabditis_elegans") ; Caenorhabditis_elegans <- Caenorhabditis_elegans # Second part is just for syntax convenience 
setName(Caenorhabditis_elegans) <- 'c_elegans'

fit_spline <- get_all_RR(all_genetic_maps = Caenorhabditis_elegans, record_names = record_names)
length(Caenorhabditis_elegans)

save(fit_spline, file = 'data/features/c_elegans_RR_spline.RData')


# 4) Melanogaster
# A specific tool will be used for D. melanogaster (RRC) and thus doesn't need a spline fit


# 5) A. thaliana (already incorporated in MareyMap)
data("Arabidopsis_thaliana") ; Arabidopsis_thaliana <-  Arabidopsis_thaliana  # Second part is just for syntax convenience 
setName(Arabidopsis_thaliana) <- 'a_thaliana'

fit_spline <- get_all_RR(all_genetic_maps = Arabidopsis_thaliana, record_names = record_names)

save(fit_spline, file = 'data/features/a_thaliana_RR_spline.RData')


# 6) S. cerevisiae
map <- read.delim("tools/RRC/s_cerevisiae/pgmap.xls", header = TRUE, sep = "\t")

# Certain markers do not have genetic distances, we thus get rid of them
to_remove <- which(is.na(map$genetic))
map <- map[-to_remove,]

split_map <- split(map, map$chr)

Saccharomyces_cerevisiae <- new("MapSet")
setName(Saccharomyces_cerevisiae) <- 'Saccharomyces_cerevisiae'

number_chromosome <- 16
# We will not use the 17th (?) chromosome, hence the 1:16
for (each_chrom in 1:number_chromosome) {
  each_map <- split_map[[each_chrom]]
  
  # The map we are using start the genetic distances at negaitve values.
  # We will shift the map so that the minimum value in the genetic distance becomes 0
  # Inspiration: https://biology.stackexchange.com/questions/35803/genetic-linkage-greater-than-50-centimorgans
  each_map$genetic <- each_map$genetic - min(each_map$genetic)
  
  # We will take the middle of the markers as physical position
  physical_positions <- unname(apply(each_map, 1, function(each_row) {
    round(mean(as.numeric(each_row[4:5])))
  })) 
  
  # Constructing the MareyMap object
  map_to_add <- new("MareyMap")
  setName(map_to_add) <- 'Saccharomyces_cerevisiae'
  mapName(map_to_add) <- as.character(each_map$chr[1])
  markerNames(map_to_add) <- as.character(each_map$gene)
  physicalPositions(map_to_add) <- physical_positions
  geneticDistances(map_to_add) <- each_map$genetic
  markerValidity(map_to_add) <- rep(TRUE, nrow(each_map))
  
  # Adding the map to the set
  Saccharomyces_cerevisiae <- Saccharomyces_cerevisiae + map_to_add
}

setName(Saccharomyces_cerevisiae) <- 's_cerevisiae'

fit_spline <- get_all_RR(all_genetic_maps = Saccharomyces_cerevisiae, record_names = record_names)

save(fit_spline, file = 'data/features/s_cerevisiae_RR_spline.RData')




