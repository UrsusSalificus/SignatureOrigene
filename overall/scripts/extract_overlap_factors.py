"""This script will extract the overall percentages of each factor, in each species
"""
from Bio import SeqIO
import sys
import glob
import os.path
import math
from joblib import Parallel, delayed

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Wanted window size:
window_size = int(sys.argv[1])
# Wanted number of threads at the same time:
n_threads = int(sys.argv[2])
# Wanted window size:
output = str(sys.argv[3])


###
# Extract all the path of files matching a certain pattern in a directory
###
def extract_path(files_directory, pattern):
    return glob.glob(files_directory + pattern)


###
# Fetch a fasta file, and clean it (remove N or n, which stands for "any nucleotides)
# Note that if the fasta file contain multiple sequences, only the first will have its CGR computed !
# Input:
#   - fasta_file : Path to the file containing the sequence one wants the CGR computed on
###
def fetch_fasta(fasta_file):
    # Will only take the first sequence of the fasta file
    try:
        records = list(SeqIO.parse(fasta_file, "fasta"))
    except:
        print("Cannot open %s, check path!" % fasta_file)
        sys.exit()
    return (records)


###
# Count all the factor only record length
###
def count_pure_record_length(record, proxies_directory):
    record_ranges = proxies_directory + '/' + record.id

    if os.path.isfile(record_ranges):
        length_factor_only = 0
        with open(record_ranges, 'r') as all_ranges:
            for each_range in all_ranges:
                line_range = each_range.strip().split()
                # Add the range length
                length_factor_only += len(range(int(line_range[0]), int(line_range[1]) + 1))
    else:
        length_factor_only = 0

    return length_factor_only


###
# Extract ranges as a list of ranges
###
def extract_ranges(ranges_file_path):
    ranges = list()
    with open(ranges_file_path, 'r') as ranges_file:
        for each_line in ranges_file:
            ranges.append([int(each_line.strip().split('\t')[0]), int(each_line.strip().split('\t')[1])])
    return ranges


###
# Find intersects between two ranges (a and b)
###
def intersects(a, b):
    overlap_ranges = []
    i = j = 0
    while i < len(a) and j < len(b):
        a_left, a_right = a[i]
        b_left, b_right = b[j]
        end_pts = sorted([a_left, a_right, b_left, b_right])
        middle = [end_pts[1], end_pts[2]]

        if a_right < b_right:
            i += 1
        else:
            j += 1

        if a_right >= b_left and b_right >= a_left:
            overlap_ranges.append(middle)

    for i in range(len(overlap_ranges) - 1):
        if overlap_ranges[i][1] in (overlap_ranges[i + 1][0], overlap_ranges[i + 1][0] - 1):
            overlap_ranges[i:i + 2] = [[overlap_ranges[i][0], overlap_ranges[i + 1][1]]]

    return overlap_ranges


###
# Extract the factors overlaps in percentage
# Inputs:
#   - genomes_directory : directory of the genomes of species
#   - main_proxies_directory : main directory of the proxies obtained through the factor_proxies script
#   - all_factors : list of all wanted factors
#   - species : name of the wanted species
# Output:
#   - A list which will contain all pairwise overlap between the factors
###
def extract_overlap_percentages(genomes_directory, main_proxies_directory, species, window_size):
    # This list will contain the percentages of genome which contain the various factors
    species_percentages = list()

    # Fetch this species records for the lengths
    records = fetch_fasta(extract_path(genomes_directory, species + '*')[0])

    # Directory of the different factor proxies of this species
    species_proxies_directory = main_proxies_directory + '/' + species

    # Extract the multiple factor directories
    all_factors = [os.path.basename(each) for each in extract_path(species_proxies_directory, '/*')]

    # We will use the whole genome as genome size
    whole_genome_length = sum([len(record) for record in records])

    for each_factor in range(len(all_factors)):
        factor_directory = species_proxies_directory + '/' + all_factors[each_factor]

        # Finally, each_factor = directory of the proxies of this species and this factor
        # Sum the length of all the pure factor sequences
        factor_record_lengths = sum([count_pure_record_length(record, factor_directory) for record in records])


        percentage = ((factor_record_lengths / whole_genome_length) * 100)

        species_percentages.append([species, os.path.basename(factor_directory), str(percentage)])

    return

# Using the config file to find which species we want to know the factor percentages
all_species = [os.path.basename(each) for each in extract_path("config/species/", '*')]

# Directory which lead to the different species proxies
main_proxies_directory = "../files/factor_proxies"

genomes_directory = "../data/genomes/"

# This list will keep all the factors percentages, in all species
all_percentages = Parallel(n_jobs=n_threads)(delayed(extract_overlap_percentages)
                                             (genomes_directory, main_proxies_directory, all_factors,
                                              all_species[each_species], window_size)
                                             for each_species in range(len(all_species)))

# The percentages may then be stored
with open(output, 'w') as outfile:
    for each_species in all_percentages:
        for each_perc in each_species:
            outfile.write('\t'.join(each_perc) + '\n')
