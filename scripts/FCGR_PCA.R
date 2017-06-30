#!/usr/bin/env Rscript
library(RColorBrewer)

### PCA analysis of a set of FCGR
# Will input the arguments:
# 1. path to the FCGRs file
# 2. window size
# 3. path to the output image (iwht its name in the path)
# 4. Plot argument: if == region -> plot with region name next ot the point
#                   else, require a path to the file containing values for region as numerics
args <- commandArgs(trailingOnly=TRUE)


FCGRS <- read.table(args[1], sep = "\t", header = TRUE)
transposed <- t(FCGRS)
R <- var(transposed)
#Eigenvectors :
U <- eigen(R)$vectors

PCA_title <- paste("PCA of", ncol(transposed) ,"FCGRs (each from", args[2], "windows)")

png(args[3], width=700, height=500, units="px")
# Prepare the text to put witht the points
# If the fourth argument = region -> will simply add the region name, and simply print the dots as they are
# Else will use the input, and plot points with size according to it.
if (args[4] == 'region') {
  text_to_add <- vector('character', ncol(transposed))
  for (i in 1:ncol(transposed)) { text_to_add[i] <- paste('region', i) }
  # Plot the result
  plot(U[,1],U[,2],pch=19,xlab="PC1",ylab="PC2",main=PCA_title)
  
  # Add text to points
  text(U[,1],U[,2], text_to_add,pos=2)
  dev.off() 
} else {
  color = read.csv(args[4], sep = '\t', header = TRUE)
  size = vector(length=ncol(color))
  for (each in 1:ncol(color)) {
    size[each] = color[,each]
  }
  text_to_add <- vector('character', ncol(transposed))
  for (i in 1:ncol(transposed)) { text_to_add[i] <- paste('region', i) }
  
  #If bigger than 10 -> must chose special brewer palette
  if (length(size) > 10) {
    # Will add more color than needed, to make sure
    colours_heat = rep(rev(brewer.pal(10, 'RdYlBu')), each=ceiling(length(size)/10))
  } else {
    colours_heat = rev(brewer.pal(size, 'RdYlBu'))
  }
  increasing_order = order(size)
  plot(U[,1],U[,2],pch=21,xlab="PC1",ylab="PC2",main=PCA_title, 
       col = 'black', bg = colours_heat[increasing_order], cex =1.5)
  #text(U[,1],U[,2], text_to_add,pos=2)
  legend_number = c(paste(max(size), '% CDS', sep=''), rep (' . ', 8) , paste(min(size), '% CDS', sep=''))
  legend(args[5],  legend_number, pch =21, pt.bg = brewer.pal(10, 'RdYlBu'), pt.cex=1.5)
  dev.off() 
}






