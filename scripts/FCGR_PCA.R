#!/usr/bin/env Rscript
library(RColorBrewer)

# TODO: clean to have generic PCA script, change the round up

### PCA analysis of a set of FCGR
# Will input the arguments:
# 1. path to the FCGRs file
# 2. window size
# 3. path to the output image (with its name in the path)
# 4. Plot argument: if == region -> plot with region name next ot the point
#                   else, require a path to the file containing values for region as numerics
args <- commandArgs(trailingOnly=TRUE)


FCGRS <- read.table(args[1], sep = "\t", header = FALSE)
# Will remove any empty column
FCGRS <- FCGRS[, colSums(is.na(FCGRS)) == 0] 
transposed <- t(FCGRS[,-1])
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
  
  # Will just plot the result ot get good x limits
  plot(U[,1],U[,2])
  x_lim = par("usr")[1:2]
  range_x <- x_lim[2] - x_lim[1]
  
  # Plot the result
  plot(U[,1],U[,2],pch=19,xlab="PC1",ylab="PC2", main=PCA_title, xlim=c(min(U[,1]) + range_x * 0.1, x_lim[2]))
  
  # Add text to points
  text(U[,1], U[,2], text_to_add, pos=2)
  dev.off() 
} else {
  feature = read.csv(args[4], sep = '\t', header = FALSE)
  
  # We will check how much we can round the result and still not lose significant differences:
  size <- feature[,2]
  good_decimal <- 1
  round_up <- TRUE
  while (round(min(size), good_decimal) == 0) {
    good_decimal = good_decimal + 1
    # If the minimum region contain less than 0.1% CDS, stop here
    if (good_decimal == 3){
      break
    }
  }
  size <- round(size*100, good_decimal)
  
  text_to_add <- vector('character', ncol(transposed))
  for (i in 1:ncol(transposed)) { text_to_add[i] <- paste('region', i) }
  
  #If bigger than 10 -> must chose special brewer palette
  if (length(size) > 10) {
    # Will add more feature than needed, to make sure
    colours_heat = rep(rev(brewer.pal(10, 'RdYlBu')), each=ceiling(length(size)/10))
  } else {
    colours_heat = rev(brewer.pal(size, 'RdYlBu'))
  }
  increasing_order = order(size)
  
  # Will just plot the result ot get good x limits
  plot(U[,1],U[,2])
  x_lim = par("usr")[1:2]
  range_x <- x_lim[2] - x_lim[1]
  
  plot(U[,1],U[,2],pch=21,xlab="PC1",ylab="PC2",main=PCA_title, xlim=c(min(U[,1]) + range_x * 0.1, x_lim[2]),
       col = 'black', bg = colours_heat[increasing_order], cex =1.5)
  #text(U[,1],U[,2], text_to_add,pos=2)
  legend_number = c(paste(max(size), '% CDS', sep=''), rep (' . ', 8) , paste(min(size), '% CDS', sep=''))
  legend('topleft',  legend_number, pch =21, pt.bg = brewer.pal(10, 'RdYlBu'), pt.cex=1.5)
  dev.off() 
}




