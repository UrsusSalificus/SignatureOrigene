#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# Quick check if we do have the directories, if not create them
function check_parent {
    parent_dir=$( dirname $1 )
    if [[ ! -d $parent_dir ]]; then
        mkdir $parent_dir
    fi
}

# Will download genetic map, and then find the fitting spline of recombination rates through a Mareymap
species="$1"
output_rr="$2"

mkdir temp
cd temp

# The downloading process will vary in between species:
if [[ $species == h_sapiens || $species == hsap_sample ]]; then
    wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
    tar -xzf genetic_map_HapMapII_GRCh37.tar.gz
    rm genetic_map_HapMapII_GRCh37.tar.gz
    cat genetic_map* >> all_genetic_maps.txt
fi

if [[ $species == m_musculus || $species == mmus_sample ]]; then
    wget http://cgd.jax.org/mousemapconverter/Revised_HSmap_SNPs.csv
fi

# Caenorhabditis elegans was already computed, as the source file is not available online.
# The genetic map of C. elegans WS185 graciously lent by MV Rockman, from his paper
# Recombinational Landscape and Population Genomics of Caenorhabditis elegans
# http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000419

if [[ $species == d_melanogaster ]]; then
    wget https://petrov.stanford.edu/RRC_scripts/RRCv2.3.tar.gz
    tar -xzf RRCv2.3.tar.gz
fi

if [[ $species == a_thaliana ]]; then
    wget https://www.arabidopsis.org/download_files/Maps/mapviewer_data/RI.data
    wget https://www.arabidopsis.org/download_files/Maps/mapviewer_data/TAIR9_AGI_marker.data
fi

# This will compute the fit spline using what we downloaded
Rscript ../scripts/MareyMap_preparation.R $species $output_rr

# Finally, we clean everything:
cd ..
rm -r temp
