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

# Will download the feature table of the species
species="$1"
# Checking if parent directory already exist
output_table="$2"
check_parent $output_table

# Using accession script to find the right accession
accession=$( bash ../input/accessions.sh $species )

# The 10th element contains the file "name", which we use to extract the files we need
genes=$(echo $accession | cut -f 10 -d '/')'_genomic.gff.gz'

wget $accession$genes
# We will add the UTR using NCBI python code
python ../scripts/add_utrs_to_gff.py $genes > temp
mv temp $output_table
# We will also add the introns using own code
python3 ../scripts/add_introns_to_gff.py $output_table temp
mv temp $output_table
# Clean everything
rm $genes