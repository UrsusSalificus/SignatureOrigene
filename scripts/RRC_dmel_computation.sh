#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# Will Calculate the Recombination Rate of a set of windows from Drosophila Melanogaster

cd ../tools/RRC/dmel/
perl RRC-open-v2.3.pl -M all_windows
mv all_windows.rrc ../../files/features/15000_RR/dmel_536.txt
rm all_windows
rm temp