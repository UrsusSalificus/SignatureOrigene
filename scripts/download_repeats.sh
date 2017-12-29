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

# Will download the RepeatMasker output of the species
species="$1"
# Checking if parent directory already exist
# (except for E. coli -> we already have the repeats for it):
if [ $species != e_coli ]; then
    output_repeats="$4"
    check_parent $output_repeats
fi

# Using accession script to find the right accession
accession=$( bash ../input/accessions.sh $species )

# Will cut the accession in array delimited by '/'
IFS='/' read -r -a array <<< $accession
# The 9th element contains the file "name", which we use to extract the files we need
if [ $species != e_coli ]; then
    repeats=$( echo ${array[9]}'_rm.out.gz' )
fi

if [ $species != e_coli ]; then
    wget $accession$repeats
    gunzip < $repeats > $output_repeats
    rm $repeats
fi
