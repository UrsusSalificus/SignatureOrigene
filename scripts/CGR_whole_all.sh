#!/usr/bin/env bash

# Name of all the different species
species='h_sapiens c_elegans d_melanogaster e_coli'
number_genomes=$( echo ${species} | wc -w )

# Find all species
for each_genome in $( seq 1 1 ${number_genomes} ) ; do
    each_species=$( echo ${species} | cut -d" " -f${each_genome} )
    species_folder=$'../data/genomes/'${each_species}

    genome_fasta_files=$( ls ${species_folder} | grep "$each_species" )
    # It's always the second one (first being cds)
    whole_genome_fasta=$( echo ${genome_fasta_files} | cut -d" " -f2 )

    infile=${species_folder}'/'${whole_genome_fasta}
    outfile=${species_folder}'/whole_CGR/each_species_CGR'

    python3 CGR_whole.py ${infile} ${outfile}
done

