#!/usr/bin/env Rscript

### PCA analysis of a set of FCGR
# Will input the args
args = commandArgs(trailingOnly=TRUE)
FCGRS <- read.table(args[1], sep = "\t", header = TRUE)
transposed <- t(FCGRS)
R <- var(transposed)
#Eigenvectors :
U <- eigen(R)$vectors

# Prepare the text to put witht the points
regions <- vector('character', ncol(transposed))
for (i in 1:ncol(transposed)) { regions[i] <- paste('region', i) }

PCA_title <- paste("PCA of", ncol(transposed) ,"FCGRs (each from", args[2], "windows)")

# Plot the result
png(args[3], width=981, height=680, units="px")
plot(U[,1],U[,2],pch=19,xlab="PC1",ylab="PC2",main=PCA_title)
text(U[,1],U[,2],regions,pos=2)
dev.off() 