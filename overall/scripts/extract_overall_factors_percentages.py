"""This script will extract the overall percentages of each factor, in each species
"""
from Bio import SeqIO
import sys
import glob
import os.path
import math

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Wanted window size:
window_size = int(sys.argv[1])
# Output file
output = str(sys.argv[2])


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


def factor_percentages(main_proxies_directory, genomes_directory, species, window_size):
    # This list will contain the percentages of genome which contain the various factors
    species_percentages = list()

    # Fetch this species records for the lengths
    records = fetch_fasta(extract_path(genomes_directory, species + '*')[0])

    # Directory of the different factor proxies of this species
    species_proxies_directory = main_proxies_directory + '/' + species

    # Extract the individual factor directories
    all_factors = [os.path.basename(each) for each in extract_path(species_proxies_directory, '/*')]

    # We will use the whole genome as genome size
    whole_genome_length = sum([len(record) for record in records])

    for each_factor in range(len(all_factors)):
        factor_directory = species_proxies_directory + '/' + all_factors[each_factor]

        # Finally, each_factor = directory of the proxies of this species and this factor
        # Sum the length of all the pure factor sequences
        factor_record_lengths = sum([count_pure_record_length(record, factor_directory) for record in records])

        # Percentage of the genome covered by this factor:
        percentage = ((factor_record_lengths / whole_genome_length) * 100)
        # If it is too low, simply count it as 0
        if percentage < 0.0001:
            percentage = 0

        # Number of windows we can get from this factor
        n_windows = int(math.floor(factor_record_lengths / window_size))

        species_percentages.append([species, os.path.basename(factor_directory), str(percentage),
                                    str(window_size), str(n_windows)])

    return species_percentages


# Using the config file to find which species we want to know the factor percentages
all_species = [os.path.basename(each) for each in extract_path("config/species/", '*')]

# Directory which lead to the different species proxies
main_proxies_directory = "../files/factor_proxies"

genomes_directory = "../data/genomes/"

# This list will keep all the factors percentages
all_percentages = list()
for each_species in all_species:
    all_percentages.extend(factor_percentages(main_proxies_directory, genomes_directory, each_species, window_size))

# The percentages may be then stored
with open(output, 'w') as outfile:
    for each_perc in all_percentages:
        outfile.write('\t'.join(each_perc) + '\n')