#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# Any argument can be passed to snakemake as argument of this setup script
# E.g. number of cores (--cores #) or dry run (-np = will just show the jobs to do, not really running them)
snakemake_arguments=$*

# Copy/paste at each step which species/features/window size should be included
# SPECIES :     h_sapiens m_musculus c_elegans d_melanogaster a_thaliana s_cerevisiae e_coli
SPECIES="h_sapiens m_musculus c_elegans d_melanogaster a_thaliana s_cerevisiae e_coli"
# WINDOW SIZE : e.g.    5000 15000 150000
WINDOWS=15000
# FEATURE TYPE :    CDS LCR
FEATURES="CDS LCR"
# KMER : e.g.       4 7
KMER="7"
#FIGURES:   e.g.    MDS correlation ratios
FIGURES="MDS correlation ratios"



# We will have to check the downloaded files, as they are input files and rise error in snakemake...
for each_species in $SPECIES; do
    for each_window in $WINDOWS; do
        for each_feature in $FEATURES; do
            for each_kmer in $KMER; do
                for each_figures in $FIGURES; do
                    if [[ $each_figures == 'ratios' ]] ; then
                        snakemake $snakemake_arguments \
                            files/results/$each_window\_$each_kmer\_ratios/$each_species.png
                    else
                        snakemake $snakemake_arguments \
                            files/results/$each_window\_$each_kmer\_$each_feature/$each_species\_$each_figures.png
                    fi
                done
            done
        done
    done
done
