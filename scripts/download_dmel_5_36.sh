#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# Will download specific version 5.36 of Drosophila Melanogaster genome from FlyBase
# The name dmel_536 was chosen to be consistent with the other species naming (something_something)
cd ../data/genomes/
wget ftp://ftp.flybase.net/releases/FB2011_04/dmel_r5.36/fasta/dmel-all-chromosome-r5.36.fasta.gz
gunzip < dmel-all-chromosome-r5.36.fasta.gz > dmel_536_genomes.fna
rm dmel-all-chromosome-r5.36.fasta.gz

# Only some records are useful for us, we will simply remove them from the fasta file
available='2L 2R 3L 3R X'
for each_available in $available; do
    # First write the header of the sequence
    echo \>"$each_available" >> temp
    # Then can fetch the sequences in  themselves using awk
    #   Awk:
    #       - /'pattern1'/ {flag=1;next}  -> Initial pattern found:
    #           -> turn on the flag and read the next line
    #       - /'pattern2'/ {flag = 0}  -> Final pattern found:
    #           -> turn off the flag
    #       - flag {print}  -> print the flagged lines
    awk '/'ID="$each_available"\;'/{flag=1;next}/>/{flag=0} flag {print}' dmel_536_genomes.fna >> temp
done
rm dmel_536_genomes.fna
mv temp dmel_536_genomes.fna

