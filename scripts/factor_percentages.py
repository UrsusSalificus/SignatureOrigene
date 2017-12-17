"""This script will extract the overall percentages of each factor
"""
from Bio import SeqIO
import sys
import glob
import os.path
import math

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Species genome path:
species_genome = str(sys.argv[1])
# Species abbreviation:
species = '_'.join(str(species_genome.split('/')[-1]).split('_')[:2])
# Wanted window size:
window_size = int(sys.argv[3])
# Output file
output = str(sys.argv[4])


###
# Check if parent directory is present, if not create it
###
def checking_parent(file_path):
    # We don't need the file name, so will take everything but the last part
    parent_directories = '/'.join(file_path.split('/')[0:(len(file_path.split('/')) - 1)])
    # As we uses parallel, we ended up we one thread doesn't seeing the directory, attempting
    # creating it, while another just did the same -> error "The file already exist", and stopped everything...
    try:
        if not os.path.exists(parent_directories):
            os.makedirs(parent_directories)
    except:
        pass


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


def factor_percentage(records, whole_genome_length, species_proxies_directory, species, factor, window_size):
    factor_directory = species_proxies_directory + '/' + factor

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

    factor_percentage = [species, os.path.basename(factor_directory), str(percentage),
                         str(window_size), str(n_windows)]

    return factor_percentage


# Fetch this species records for the lengths
records = fetch_fasta(species_genome)

# We will use the whole genome as genome size
whole_genome_length = sum([len(record) for record in records])

# Directory of the different factor proxies of this species
species_proxies_directory = "../files/factor_proxies" + '/' + species

# Extract the individual factor directories
all_factors = [os.path.basename(each) for each in extract_path(species_proxies_directory, '/*')]

# This list will contain the percentages of genome which contain the various factors
species_percentages = list()

for each_factor in range(len(all_factors)):
    species_percentages.append(factor_percentage(records, whole_genome_length, species_proxies_directory, species,
                                                 all_factors[each_factor], window_size))

checking_parent(output)
# The percentages may be then stored
with open(output, 'w') as outfile:
    for each_perc in species_percentages:
        outfile.write('\t'.join(each_perc) + '\n')
