#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# This script concatenate into one file all the pure FCGRs of the wanted list of factors in two different species
windows="$1"
n_samples="$2"
kmer="$3"
species="$4"
all_factors="$5"
output="$6"


function get_factor {
    file_name=$( basename "$1" )
    IFS='_' read -r -a array <<< $file_name
    echo ${array[0]}
}

# Extract all the factors
sep_factors=$( echo $all_factors | tr _ ' ' )
# Add uncategorized
sep_factors+=" uncategorized"

# Find all pure factors FCGRS of this species
all_species_files=$( find files/FCGRs/$windows\_$n_samples\_$kmer/$species/*_pure_FCGRs.txt )

# Extract only the factors we want
for each_file in $all_species_files; do
    for each_factor in $sep_factors; do
        if [[ "$( get_factor $each_file )" == "$each_factor" ]]; then
            files_kept+="$each_file "
        fi
    done
done

# Now, can concatenate everything in one file
cat $files_kept ../scaling/files/FCGRs/$windows\_$n_samples\_$kmer/$species\_FCGRs.txt >> $output
