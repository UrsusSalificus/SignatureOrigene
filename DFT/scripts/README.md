### Scripts used for the analysis
#### Not up to date...

Brief explanation of files found in this directory : 

* `all_windowed_CGR.py`  :  for each species, for each window of a certain size, compute the 
Chaos Game Representation (CGR) of a sequence.
* `all_windowed_FCGR.py`  :  For each CGR of a certain window size, compute all the different k-mer **F**requencies 
* `CGR_functions.py`  :  gathers all the different functions used to compute and analyse the Chaos Game Representation 
(CGR) of a sequence
* `CGR_plot.R` : going from the coordinates to the actual CGR picture 
* `download_genomes.sh` : download all the genomes used for the analysis, and untar them
* `FCGR_PCA.R` : output a picture of PC1 VS PC2 of a set of FCGR
* `impact_of_N.py` : overview of the impact of either removing (which creates "wrong" k-mer) unknown nucleotides (N),
or ignoring any window containing an unknown nucleotide
*  `bsub*` : LSF jobs of the different python scripts




