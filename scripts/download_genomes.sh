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

# Will download whole concatenated genomes
species="$1"
output_genome="$2" ; check_parent $output_genome
output_table="$3" ; check_parent $output_table
# We already have the repeats data of E. coli:
if [ $species != e_coli ]; then
    output_repeats="$4" ; check_parent $output_repeats
fi

# Using hash table to find hte right accession path
declare -A all_accessions
all_accessions=(["h_sapiens"]="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.36_GRCh38.p10/"
["hsap_sample"]="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.36_GRCh38.p10/"
["m_musculus"]="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.25_GRCm38.p5/"
["mmus_sample"]="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.25_GRCm38.p5/"
["c_elegans"]="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/"
["d_melanogaster"]="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/"
["a_thaliana"]="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.3_TAIR10/"
["s_cerevisiae"]="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/"
["e_coli"]="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/")

accession=${all_accessions[$species]}

# Will cut the accession in array delimited by '/'
IFS='/' read -r -a array <<< $accession
# The 9th element contains the file "name", which we use to extract the files we need
whole_genome=$( echo ${array[9]}'_genomic.fna.gz' )
feature_table=$( echo ${array[9]}'_feature_table.txt.gz' )

if [ $species != e_coli ]; then
    repeats=$( echo ${array[9]}'_rm.out.gz' )
fi

wget $accession$whole_genome
gunzip < $whole_genome > $output_genome
rm $whole_genome

wget $accession$feature_table
gunzip < $feature_table > $output_table
rm $feature_table

if [ $species != e_coli ]; then
    wget $accession$repeats
    gunzip < $repeats > $output_repeats
    rm $repeats
fi

if [ $species == "hsap_sample" ] || [ $species == "mmus_sample" ]; then
    python3 scripts/sampling.py $output_genome
fi