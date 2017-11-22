#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# Determine which python script to use to extract the wanted feature

factor="$1"
window_size="$2"
repeats="$5"
features="$6"
species="$9"
all_arguments="$*"

if [[ $factor == 'RR' ]]; then
    python3 scripts/extract_RR.py "$all_arguments"
else
    # build follow_up path
    follow_up="../data/following/factor_proxies/$window_size/$species\_$factor\_proxies_done.txt"
    # Extract the proxy ranges
    python3 ../scripts/factor_proxies.py "$1" "$2" "$3" "$5" "$6" "$follow_up"
    # Extract the percentages using the ranges
    python3 scripts/extracting_factor_percentage.py "$all_arguments"
fi