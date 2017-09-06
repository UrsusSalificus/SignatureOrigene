#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

library(data.table)
library(ggplot2)

### Nucleotides ratios analyis
# Will input the arguments:
# 1. path to output
# 2. path to the distance matrix
# 3. path to the ratio table
# 4. feature type
# 5. window size
# 6. Genomic signature type
# 7. IF USING FCGR : k-mer size

output <- args[1]
gs <- basename(args[6])
feature_type <- args[4]
window_size <- args[5]
if (gs == 'FCGRs'){
  kmer <- args[7]
} 

distance_matrix <- readRDS(args[2])
distance_matrix <- as.matrix(distance_matrix)

feature <- fread(args[3], sep = "\t", header = TRUE)
# Remove last column (empty column)
feature <- feature[,1:(ncol(feature)-1)]

mean_dist <- rowMeans(distance_matrix)
far_away <- order(mean_dist)[round(length(mean_dist)*0.95):length(mean_dist)]
far <- rep('FALSE', time = length(mean_dist))
far[far_away] <- 'TRUE'

feature$far <- far

feature = melt(feature, id.vars = c("record", "far"),
               measure.vars = c("A", "C", "T", "G", "AG", "CG"))

ggplot(data = feature, aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill = far), width = 1) + theme_bw()

by(feature, feature$variable, function (each_ratio) {
  wilcox.test(each_ratio$value[each_ratio$far == "TRUE"], each_ratio$value[each_ratio$far == "FALSE"], paired=FALSE)
})

test <- feature[feature$variable == 'A']
lol <- order(test$value)
ggplot(data = test, aes(x = c(1:309), y = value[lol])) + 
  geom_point(aes(col = far[lol])) + theme_bw()



