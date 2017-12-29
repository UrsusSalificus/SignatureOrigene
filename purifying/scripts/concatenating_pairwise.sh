#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# This script concatenate into one file all the pure FCGRs of the wanted list of factors in two different species
windows="$1"
n_samples="$2"
kmer="$3"
species="$4"
comparison="$5"
all_factors="$6"
output="$7"


function get_factor {
    file_name=$( basename "$1" )
    IFS='_' read -r -a array <<< $file_name
    echo ${array[0]}
}


# This directory will contain all our files temporarily
mkdir temp

# Extract all the factors
sep_factors=$( echo $all_factors | tr _ ' ' )

# Find all pure factors FCGRS of this species
all_species_files=$( find files/FCGRs/$windows\_$n_samples\_$kmer/$species/*_pure_FCGRs.txt )

# Extract only the factors we want
for each_file in $all_species_files; do
    for each_factor in $sep_factors; do
        if [[ "$( get_factor $each_file )" == "$each_factor" ]]; then
            all_actual_species+="lol_$each_factor "
        fi
    done
done

# Now store the right files in temp
for each_actual in $all_actual_species; do
    file_name_only=$( basename $each_actual )
    # Add the species name at the start of each line
    sed -e "s/^/$species\_/" $each_actual > temp/$species\_$file_name_only
done


# Same thing for the comparison
all_comparison_files=$( find files/FCGRs/$windows\_$n_samples\_$kmer/$comparison/*_pure_FCGRs.txt )

# Extract only the factors we want
for each_file in $all_comparison_files; do
    for each_factor in $sep_factors; do
        if [[ "$( get_factor $each_file )" == "$each_factor" ]]; then
            all_comparison+="lol_$each_factor "
        fi
    done
done

for each_comparison in $all_comparison; do
    file_name_only=$( basename $each_comparison )
    # Add the species name at the start of each line
    sed -e "s/^/$comparison\_/" $each_comparison > temp/$comparison\_$file_name_only
done


# Now, can concatenate everything in one file
cat temp/$species*.txt temp/$comparison*.txt >> $output
rm -r temp