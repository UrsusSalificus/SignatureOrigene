#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# Any argument can be passed to snakemake as argument of this setup script
# E.g. number of cores (--cores #) or dry run (-np = will just show the jobs to do, not really running them)
snakemake_arguments=$*

# Copy/paste at each step which species/features/window size should be included
# SPECIES :    h_sapiens hsap_sample m_musculus mmus_sample c_elegans d_melanogaster a_thaliana s_cerevisiae e_coli
SPECIES="hsap_sample mmus_sample c_elegans d_melanogaster a_thaliana s_cerevisiae e_coli"
# WINDOW SIZE : e.g.    5000 15000 150000
WINDOWS=15000
# FEATURE TYPE :    CDS RNA LCR TE tandem
FEATURES="CDS RNA LCR TE tandem"
# KMER : e.g.       4 7
KMER="7"


for each_species in $SPECIES; do
  # Launching the whole Snakemake cascade
    for each_window in $WINDOWS; do
        for each_kmer in $KMER; do
            # A) We need the whole genome distance matrix, which will be computed through the scaling snakemake
            cd ../scaling
            snakemake $snakemake_arguments \
                files/distances/pearson/$each_window\_$each_kmer/$each_species\_dist_matrix.RData
            cd $go_back

            # B) We now will compute both the masked and pure factors distance matrix VS center
            for each_factor in $FACTORS; do
                # B.1) Masked will be done all on the masking directory
                snakemake $snakemake_arguments \
                    files/distances/pearson/$each_window\_$each_kmer/$each_species\_$each_factor\_masked_vs_center_dist_matrix.RData

                # B.2) Pure will need the FCGRs found in the purifying directory
                cd ../purifying
                snakemake $snakemake_arguments \
                    files/FCGRs/$each_window\_$each_kmer/$each_species\_$each_factor\_pure_FCGRs.txt
                cd $go_back
                # Which will then be used to find the pure VS center distance matrix
                snakemake $snakemake_arguments \
                    files/distances/pearson/$each_window\_$each_kmer/$each_species\_$each_factor\_pure_vs_center_dist_matrix.RData
            done

            # C) Our control will be the whole genome VS center distance
            snakemake $snakemake_arguments \
                    files/distances/pearson/$each_window\_$each_kmer/$each_species\_whole_vs_center_dist_matrix.RData
            # TODO: find why we have to recompute the distance this way instead of using the whole distance matrix

            # D) We would finally be able to use both distance matrices to produces boxplots
            snakemake $snakemake_arguments \
                    files/results/$each_window\_$each_kmer/$each_species\_boxplots_all_factors.png
        done
    done
done
