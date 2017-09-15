#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# Will download chromosome 8 of specific version GRCh36 of Homo sapiens genome
# The name hsap_36sample was chosen to be consistent with the other species naming (something_something)
cd ../data/genomes/
wget http://ftp.ensembl.org/pub/release-54/fasta/homo_sapiens/dna/Homo_sapiens.NCBI36.54.dna.chromosome.10.fa.gz
gunzip < Homo_sapiens.NCBI36.54.dna.chromosome.10.fa.gz > hsap_36sample_genomes.fna
rm Homo_sapiens.NCBI36.54.dna.chromosome.10.fa.gz

# Will also download the Recombination rate from HapMap release 22 (obtained using GRCh36)
cd ../../tools/RRC/hsap/
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/latest/rates/genetic_map_chr10_b36.txt

### TO REMOVE
#snakemake --cores 3 files/distances/pearson/15000_7/hsap_36sample_dist_matrix.txt
