#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2018 Titouan Laessle
# License : MIT

# Will reproduce the whole analysis performed for the master thesis
# Results can be found in purifying/files/results and masking/files/results directories

# Any argument can be passed to snakemake as argument of this script
# E.g. number of cores (--cores #) or dry run (-np = will just show the jobs to do, not really running them)
snakemake_arguments="$*"

# Set the different inputs
each_window=15000
each_sample=100
each_kmer=7
SPECIES="h_sapiens m_musculus c_elegans d_melanogaster a_thaliana s_cerevisiae e_coli"
FACTORS="CDS intron UTR RNA LCR TE tandem"
all_factors="CDS_intron_UTR_RNA_LCR_TE_tandem"


### First part: MDS of pure features ###

# Specific inputs
each_figure="MDS"

cd purifying
# Remember where we are :
go_back=$( pwd )

# The first part is to extract and clean factor ranges:
for each_species in $SPECIES; do
    for each_factor in $FACTORS; do
        # If feature = UTR, do not compute for S. cerevisiae and E. coli
        if [[ $each_factor == 'UTR' ]] && \
            ([[ $each_species == 's_cerevisiae' ]] || [[ $each_species == 'e_coli' ]]); then
            echo "$each_factor not computed for $each_species"
        elif [[ $each_factor == 'intron' && $each_species == 'e_coli' ]] ; then
            echo "$each_factor not computed for $each_species"
        else
            # Finding factor ranges
            snakemake $snakemake_arguments \
                ../data/following/factor_proxies/$each_species/$each_factor\_proxies_done.txt
        fi
    done
    # After finding all the factors' ranges, we must clean them from overlaps
    snakemake $snakemake_arguments \
        ../data/following/factor_filtered/$each_window/$each_species\_done.txt
done

# The second part goes from these ranges to find sequences where we only find the factor
for each_species in $SPECIES; do
    # We now have new filtered factors
    cd ../config/new_factors/$each_species
    NEW_FACTORS=$( find * )
    cd $go_back

    # A) This part will compute all the FCGRs of pure sequences we need from the purifying snakemake
    for each_factor in $NEW_FACTORS; do
        # If feature = UTR, do not compute for S. cerevisiae and E. coli
        if [[ $each_factor == 'UTR' ]] && \
            ([[ $each_species == 's_cerevisiae' ]] || [[ $each_species == 'e_coli' ]]); then
            echo "$each_factor not computed for $each_species"
        elif [[ $each_factor == 'intron' && $each_species == 'e_coli' ]] ; then
            echo "$each_factor not computed for $each_species"
        else
            snakemake $snakemake_arguments \
                files/FCGRs/$each_window\_$each_sample\_$each_kmer/$each_species/$each_factor\_pure_FCGRs.txt
        fi
    done

    # B) We also need the whole genome FCGRs, which will be computed through the scaling snakemake
    go_back=$( pwd )
    cd ../scaling
    snakemake $snakemake_arguments \
        files/FCGRs/$each_window\_$each_sample\_$each_kmer/$each_species\_FCGRs.txt
    cd $go_back
done

# C) Figures
for each_species in $SPECIES; do
    # From the concatenated FCGRs of all factors in each species, to the MDS these factors
    # (concatenating FCGRs of factor + whole -> distance matrix -> fitting -> MDS)
    snakemake $snakemake_arguments \
        files/results/$each_window\_$each_sample\_$each_kmer/$all_factors/$each_species\_MDS_all_factors.png
done

cd $go_back
cd ..



### Second part: comparison of distance to center ###

cd masking
# Remember where we are :
go_back=$( pwd )

# The second part uses ranges previously computed to find sequences where the factor is absent
for each_species in $SPECIES; do
    # A) We need the whole genome distance matrix, which will be computed through the scaling snakemake
    cd ../scaling
    snakemake $snakemake_arguments \
        files/distances/manhattan/$each_window\_$each_sample\_$each_kmer/$each_species\_dist_matrix.RData
    cd $go_back

    # B) We now will compute both the masked and pure factors distance matrix VS center
    for each_factor in $NEW_FACTORS; do
        # If feature = UTR, do not compute for S. cerevisiae and E. coli
        if [[ $each_factor == 'UTR' ]] && \
            ([[ $each_species == 's_cerevisiae' ]] || [[ $each_species == 'e_coli' ]]); then
            echo "$each_factor not computed for $each_species"
        elif [[ $each_factor == 'intron' && $each_species == 'e_coli' ]] ; then
            echo "$each_factor not computed for $each_species"
        else
            # B.1) Masked will be done all on the masking directory
            snakemake $snakemake_arguments \
                files/distances/manhattan/$each_window\_$each_sample\_$each_kmer/$each_species/$each_factor\_masked_vs_center_dist_matrix.RData

            # B.2) Pure will need the FCGRs found in the purifying directory
            cd ../purifying
            snakemake $snakemake_arguments \
                files/FCGRs/$each_window\_$each_sample\_$each_kmer/$each_species/$each_factor\_pure_FCGRs.txt
            cd $go_back
            # B.3) Which will then be used to find the pure VS center distance matrix
            snakemake $snakemake_arguments \
                files/distances/manhattan/$each_window\_$each_sample\_$each_kmer/$each_species/$each_factor\_pure_vs_center_dist_matrix.RData
        fi
    done

    # C) Our control will be the whole genome VS center distance
    snakemake $snakemake_arguments \
            files/distances/manhattan/$each_window\_$each_sample\_$each_kmer/$each_species/whole_vs_center_dist_matrix.RData
    # TODO: find why we have to recompute the distance this way instead of using the whole distance matrix

    # D) We would finally be able to use both distance matrices to produces boxplots
    snakemake $snakemake_arguments \
            files/results/$each_window\_$each_sample\_$each_kmer/$all_factors/$each_species\_boxplots_all_factors.png
done

cd $go_back
cd ..



### Third part: the human/mouse comparison ###

SPECIES="h_sapiens m_musculus"

cd purifying

# Basically, all the FCGRs are computed, now need to concatenate these factors ($all_factors) in the two species
for each_species in $SPECIES; do
    # Comparing the concatenated FCGRs of all factors between species
    # Only fo this if it doesn't already exist (reverse order of species/comparison)
    for each_comparison in $SPECIES; do
        if [[ $each_comparison != $each_species ]] && \
            [[ ! -f files/results/$each_window\_$each_sample\_$each_kmer/$all_factors/$each_comparison\_vs_$each_species\_MDS_pairwise.png ]] ; then
            snakemake $snakemake_arguments \
            files/results/$each_window\_$each_sample\_$each_kmer/$all_factors/$each_species\_vs_$each_comparison\_MDS_pairwise.png
        fi
    done
done



### Fourth part: the human/mouse comparison but without tandem and LCR ###

all_factors="CDS_intron_UTR_RNA_TE"

# Again, all the FCGRs are already computed, the difference resides in what to concatenates ($all_factors)
for each_species in $SPECIES; do
    # Comparing the concatenated FCGRs of all factors between species
    # Only fo this if it doesn't already exist (reverse order of species/comparison)
    for each_comparison in $SPECIES; do
        if [[ $each_comparison != $each_species ]] && \
            [[ ! -f files/results/$each_window\_$each_sample\_$each_kmer/$all_factors/$each_comparison\_vs_$each_species\_MDS_pairwise.png ]] ; then
            snakemake $snakemake_arguments \
            files/results/$each_window\_$each_sample\_$each_kmer/$all_factors/$each_species\_vs_$each_comparison\_MDS_pairwise.png
        fi
    done
done



### Again, results can now be found in purifying/files/results and masking/files/results directories ###