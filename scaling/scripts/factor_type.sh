#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# Determine which python script ot use to extract the wanted feature

factor="$1"
other_arguments="$*"

if [[ $factor == 'RR' ]]; then
    python3 scripts/extract_RR.py $other_arguments
else
    python3 scripts/extracting_factor_percentage.py $other_arguments
fi