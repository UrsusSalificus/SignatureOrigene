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
# 4. window size
# 5. Genomic signature type
# 6. IF USING FCGR : k-mer size
args <- commandArgs(trailingOnly=TRUE)

output <- args[1]
window_size <- args[4]
gs <- basename(args[5])
if (gs == 'FCGRs'){
  kmer <- args[6]
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
               measure.vars = c("A", "C", "T", "G", "AG", "CG", "TG"))

# To enable space at the top of the boxplots to put if significant
py <- pretty(feature$value)

if (gs == 'FCGRs'){
  plot_title <- paste('Comparison of ratios of nucleotides\n',
                      'Using ', gs, ' with k = ', kmer, ' and ', window_size, ' bp windows ', sep = '')
} else {
  plot_title <- paste('Comparison of ratios of nucleotides\n',
                      'Using ', gs, ' with ', window_size, ' bp windows', sep = '')
}


g <- ggplot(data = feature, aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill = far), width = 1) + theme_bw() +
  scale_y_continuous(breaks=py, limits=range(py)) + 
  scale_fill_discrete(name="Distance to\nother windows",
                      labels=c("Central", "Peripherical")) + 
  labs(title=plot_title, y = "Ratios [count/total]") +
  theme(
    plot.title = element_text(size = 15, face="bold"),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 12),
    axis.title.y = element_text(size = 15),
    axis.text.y  = element_text(size = 12),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12)
  )

p_val <- as.numeric(by(feature, feature$variable, function (each_ratio) {
  wilcox.test(each_ratio$value[each_ratio$far == "TRUE"], each_ratio$value[each_ratio$far == "FALSE"], paired=FALSE)$p.value
}))
p_val <- p.adjust(p_val, method = 'holm')
significant <- which(p_val <= 0.05)

good_spots <- sapply(significant, function (each_sign) {
  each_ratio <- levels(as.factor(feature$variable))[each_sign]
  range_ratio <- max(feature$value) - min(feature$value)
  max_ratio <- max(feature$value[feature$variable == each_ratio])
  good_spot <- max_ratio + (range_ratio * 0.1)
  return(good_spot)
})

code <- sapply(significant, function (each_sign){
  if (p_val[each_sign] < 0.001) {
    return('***')
  } else if (p_val[each_sign] < 0.01) {
    return('**')
  } else {
    return('*')
  }
})


g <- g + annotate("text", x = significant, y = good_spots, label = code, size = 8)

png(output, width=700, height=500, units="px")
print(g)
dev.off() 



