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
feature_table=$( echo $accession | cut -f 10 -d '/')'_feature_table.txt.gz'

wget $accession$feature_table
gunzip < $feature_table > $output_table
rm $feature_table