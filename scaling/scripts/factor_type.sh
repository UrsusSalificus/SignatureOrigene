#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# Determine which python script to use to extract the wanted feature

factor="$1"
window_size="$2"
species_genome="$3"
IFS='_' read -r -a array <<< $( basename "$species_genome" )
species=$( echo ${array[0]}_${array[1]} )
species_sample="$4"
repeats="$5"
features="$6"
genes="$7"
RR="$8"
output="$9"

if [[ $factor == 'RR' ]]; then
    python3 scripts/extract_RR.py "$window_size" "$species_sample" "$RR" "$output"
else
    # build follow_up path
    follow_up="../data/following/factor_proxies/$window_size/$species\_$factor\_proxies_done.txt"
    # Extract the proxy ranges
    python3 ../scripts/factor_proxies.py "$factor" "$window_size" "$species_genome" "$repeats" "$features" \
        "$genes" "$follow_up"
    # Extract the percentages using the ranges
    python3 scripts/extracting_factor_percentage.py "$factor" "$window_size" "$species_genome" \
        "$species_sample" "$output"
fi