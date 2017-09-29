#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(MareyMap)
library(data.table)

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
RR_directory <- 'data/features/recombination_rates'
dir.create(RR_directory)


###################################################################################################


# 1) H. sapiens
map <- fread('tools/RRC/h_sapiens/all_genetic_maps.txt', header = TRUE)

# Will remove the unwanted headers and also remove small parts of X
map <- map[! (map$Chromosome == 'Chromosome' | map$Chromosome == 'chrX_par1' | map$Chromosome == 'chrX_par2'),]

# We will also simplify the chromosome column (from character e.g. chr1 to numeric 1)
map$Chromosome <- as.factor(map$Chromosome)
levels(map$Chromosome) <- c(1, 10:19, 2, 20:22, 3:9, 23)
map$Chromosome <- as.numeric(map$Chromosome)

split_map <- split(map, map$Chromosome)

fit_spline <- lapply(split_map, function(each_map) {
  spline <- smooth.spline(each_map$`Position(bp)`, each_map$`Rate(cM/Mb)`, cv= TRUE)
  return(spline)
})

# Add the ids
chromosome_ids <- record_names$V1[record_names$V2 == 'h_sapiens']

# Associate the id to each spline fit
names(fit_spline) <- chromosome_ids

save(fit_spline, file = paste(RR_directory, 'h_sapiens_RR_spline.RData', sep = '/'))


###################################################################################################


# 2) M. musculus
map <- read.csv('tools/RRC/m_musculus/Revised_HSmap_SNPs.csv')

split_map <- split(map, map$chr)

Mus_musculus <- new("MapSet")
setName(Mus_musculus) <- 'Mus_musculus'

number_chromosome <- length(split_map)
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

save(fit_spline, file = paste(RR_directory, 'm_musculus_RR_spline.RData', sep = '/'))


###################################################################################################


# 3) C. elegans (already incorporated in MareyMap) 
map <- read.csv('tools/RRC/c_elegans/CelegansRIAILmap.csv')

# We will need dummy marker names
markers <- paste('marker', 1:nrow(map), sep = '')
map$markers <- markers

split_map <- split(map, map$chr)

Caenorhabditis_elegans <- new("MapSet")
setName(Caenorhabditis_elegans) <- 'Caenorhabditis_elegans'

number_chromosome <- length(split_map)
for (each_chrom in 1:number_chromosome) {
  each_map <- split_map[[each_chrom]]
  
  # Constructing the MareyMap object
  map_to_add <- new("MareyMap")
  setName(map_to_add) <- 'Caenorhabditis_elegans'
  markerNames(map_to_add) <- as.character(each_map$markers)
  physicalPositions(map_to_add) <- as.numeric(each_map$bps)
  geneticDistances(map_to_add) <- as.numeric(each_map$pos)
  mapName(map_to_add) <- as.character(each_map$chr[1])
  markerValidity(map_to_add) <- rep(TRUE, nrow(each_map))
  
  # Adding the map to the set
  Caenorhabditis_elegans <- Caenorhabditis_elegans + map_to_add
}

setName(Caenorhabditis_elegans) <- 'c_elegans'

fit_spline <- get_all_RR(all_genetic_maps = Caenorhabditis_elegans, record_names = record_names)

save(fit_spline, file = paste(RR_directory, 'c_elegans_RR_spline.RData', sep = '/'))


###################################################################################################


# 4) Melanogaster
chrX <- read.table('tools/RRC/d_melanogaster/Comeron_tables/Comeron_100kb_chrX.txt')
chr2L <- read.table('tools/RRC/d_melanogaster/Comeron_tables/Comeron_100kb_chr2L.txt')
chr2R <- read.table('tools/RRC/d_melanogaster/Comeron_tables/Comeron_100kb_chr2R.txt')
chr3L <- read.table('tools/RRC/d_melanogaster/Comeron_tables/Comeron_100kb_chr3L.txt')
chr3R <- read.table('tools/RRC/d_melanogaster/Comeron_tables/Comeron_100kb_chr3R.txt')

# Following the wanted_records.txt order
split_map <- list(chrX, chr2L, chr2R, chr3L, chr3R)

fit_spline <- lapply(split_map, function(each_map) {
  spline <- smooth.spline(each_map$V1, each_map$V2, cv= TRUE)
  return(spline)
})

# Add the ids
chromosome_ids <- record_names$V1[record_names$V2 == 'd_melanogaster']

# Associate the id to each spline fit
names(fit_spline) <- chromosome_ids

save(fit_spline, file = paste(RR_directory, 'd_melanogaster_RR_spline.RData', sep = '/'))


###################################################################################################


# 5) A. thaliana
gen_distances <- read.table('tools/RRC/a_thaliana/RI.data') 
phys_pos <- read.table('tools/RRC/a_thaliana/TAIR9_AGI_marker.data') 

# We will clean any dupplicates in the data
phys_pos <- phys_pos[!duplicated(phys_pos$V1),]

gen_distances$V1 <- as.character(gen_distances$V1)
phys_pos$V1 <- as.character(phys_pos$V1)

# Due to type conversion issues, we had to use a for loop (forgive me, Oh Lord of's Apply)
keeping <- list()
average_position <- list()
for (each_row in 1:nrow(gen_distances)){
  gen_row <- gen_distances[each_row, ]
  if (any(phys_pos$V1 == as.character(gen_row[1]))) {
    keeping[[each_row]] <- TRUE
    corresponding_row <- phys_pos[phys_pos$V1 == as.character(gen_row[1]),]
    average_position[[each_row]] <- round(mean(c(as.numeric(corresponding_row[3]), as.numeric(corresponding_row[4]))))
  } else {
    keeping[[each_row]] <- FALSE
    average_position[[each_row]] <- 0
  }
}

gen_distances$keeping <- unlist(keeping) ; gen_distances$average_pos <- unlist(average_position)

map <- gen_distances[gen_distances$keeping,]
names(map)[1:4] <- c('id', 'marker', 'genetic', 'chr')

split_map <- split(map, map$chr)

Arabidopsis_thaliana <- new("MapSet")
setName(Arabidopsis_thaliana) <- 'Arabidopsis_thaliana'

number_chromosome <- length(split_map)
for (each_chrom in 1:number_chromosome) {
  each_map <- split_map[[each_chrom]]
  
  # Sort it by genetic distance (cM):
  each_map <- each_map[order(each_map$genetic),]
  
  # Constructing the MareyMap object
  map_to_add <- new("MareyMap")
  setName(map_to_add) <- 'Arabidopsis_thaliana'
  mapName(map_to_add) <- as.character(each_map$chr[1])
  markerNames(map_to_add) <- as.character(each_map$marker)
  physicalPositions(map_to_add) <- each_map$average_pos
  geneticDistances(map_to_add) <- each_map$genetic
  markerValidity(map_to_add) <- rep(TRUE, nrow(each_map))
  
  # Adding the map to the set
  Arabidopsis_thaliana <- Arabidopsis_thaliana + map_to_add
}

setName(Arabidopsis_thaliana) <- 'a_thaliana'

fit_spline <- get_all_RR(all_genetic_maps = Arabidopsis_thaliana, record_names = record_names)

save(fit_spline, file = paste(RR_directory, 'a_thaliana_RR_spline.RData', sep = '/'))


###################################################################################################


# 6) S. cerevisiae
map <- read.delim("tools/RRC/s_cerevisiae/pgmap.xls", header = TRUE, sep = "\t")

# Certain markers do not have genetic distances, we thus get rid of them
to_remove <- which(is.na(map$genetic))
map <- map[-to_remove,]

split_map <- split(map, map$chr)

Saccharomyces_cerevisiae <- new("MapSet")
setName(Saccharomyces_cerevisiae) <- 'Saccharomyces_cerevisiae'

number_chromosome <- length(split_map)
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

save(fit_spline, file = paste(RR_directory, 's_cerevisiae_RR_spline.RData', sep = '/'))




