##  Trying to distinguish properties and origin of Genomic Signatures

This project is the result of my Master project in MLS Bioinformatic at Lausanne University (Swiss), in the  
[group of Marc Robinson-Rechavi](https://www.unil.ch/dee/robinson-rechavi-group), 
and under the supervision of [Kamil Jaron](https://github.com/KamilSJaron).

To perform the whole analysis as it was performed for the master thesis, one can simply run the 
`master_analysis.sh` script on the root of this repository.  
&#9888; **THIS MAY TAKE SEVERAL HOURS TO COMPLETE** &#9888;  
To speed up the analysis multiple cores may be allocated, using the following syntax (for 2 cores):  
`bash master_analysis.sh --cores 2`

### Description of the project
The project is very briefly presented in the wiki page [The project](https://github.com/UrsusSalificus/SignatureOrigene/wiki/The-project) 
The rest of the [wiki](https://github.com/UrsusSalificus/SignatureOrigene/wiki) also contain some useful information
such as the limitations, or how to add a new species to the analysis.

 
### I want to redo it all!
The project is built with multiple `Snakefiles`, files which are basically a big recipe, building from scratch all 
the ingredients it needs to produce a final output, in our case figures. 
But as we usually compare more than one species,
the snakefiles are themselves called by a bash script (`setup.sh`) in a for loop. 
These script are themselves called by the user, and offers many way to customise the analysis.

**Note**: to give Snakemake arguments to the the Snakefiles, one must give them as argument of `setup.sh`!

The analysis is separated in various directories, which basically treat the inputs (the genomes), differently:
- `purifying`: this directory will output a Multidimensional scaling (MDS) of the distance separating the k-mer 
frequencies of sequences only composed of features. In other words, it is a way to compare how (di)similar are 
the k-mer frequencies of each features.
- `masking`: this directory will use the k-mer frequencies of both the sequences obtained in `purifying` and
sequences which are specifically built to contain no nucleotides of the feature. 
These k-mer frequencies are then compared to a completely subjective representative k-mer frequencies of the genome.
The output is boxplots of the distances in between pure/masked sequences and this representative k-mer frequencies.
- `overall`:  this directory will output the percentages of genome annotated with each feature.
- `random`: this directory will very roughly compare the k-mer frequencies of random samples from various species.
- `data`: will be mainly built during analysis, but already contain some necessary data.
- `input`: contain various files used during configuration/plotting.

Directories contain their specific scripts, yet some general scripts are found in the `scripts` directory.

Some directories are not used anymore, and may not work properly at the moment:
- `scaling`: initially used to compare the percentage of feature found in certain windows to how similar are 
the 7-mer frequencies of each window to the representative k-mer frequencies.
- `DFT`: initially used to build Power Spectrum of windows (instead of k-mer frequencies).
  
**tl;dr**

I want the colourful dotty plot, and I want to use 2 cores whenever possible!
```
cd purifying
bash setup.sh --cores 2
```

### Dependencies
The analysis is based on [Snakemake](http://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Python 3 modules used during the analyses:
```
Bio (Biopython)
joblib
re
glob
subprocess
sys
os
math
numpy
itertools
random
```

R modules used during the analysis
```
ggplot2
MareyMap
data.table
dunn.test
igraph
grid
ggrepel
RColorBrewer
wesanderson
```
