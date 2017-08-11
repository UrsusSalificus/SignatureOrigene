#!/usr/bin/env Rscript

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

### Plotting CGR
# Stock the arguments given 
# 1. The path to the CGR file
# 2. The path of the output file (with its name)
args <- commandArgs(trailingOnly=TRUE)
# Import the CGR
dat <- read.table(args[1], sep='\t')
# Plot without any border/axis
png(args[2], width=500, height=500, units="px")
plot(dat$V1, dat$V2, pch='.', xaxt='n', yaxt = 'n', ylab='', xlab='', bty='n')
# Add the nucleotides at the corners
par(xpd=TRUE)
text(-0.05, -0.05, "A", cex = 3)
text(-0.05, 1.05, "C", cex = 3)
text(1.05, -0.05, "T", cex = 3)
text(1.05, 1.05, "G", cex = 3)
par(xpd=FALSE) 
dev.off() 


