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

# Will download whole concatenated genome of the species
species="$1"
# Checking if parent directory already exist
output_genome="$2"
check_parent $output_genome

# Using accession script to find the right accession
accession=$( bash ../input/accessions.sh $species )

# The 10th element contains the file "name", which we use to extract the files we need
whole_genome=$(echo $accession | cut -f 10 -d '/')'_genomic.fna.gz'

wget $accession$whole_genome
gunzip < $whole_genome > $output_genome
rm $whole_genome

if [[ $species == hsap_sample || $species == mmus_sample ]]; then
    python3 ../scripts/sampling_chromosomes.py $output_genome
fi

# Keep only wanted records
python3 ../scripts/keep_wanted_records.py $output_genome
