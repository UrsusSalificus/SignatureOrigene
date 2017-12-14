#!/usr/bin/env python3

"""This script will extract the uncategorized nucleotides
"""
from Bio import SeqIO
import sys
import glob
import os.path
from joblib import Parallel, delayed
import itertools


__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Species genome path:
species_genome = str(sys.argv[1])
# Species abbreviation:
species = '_'.join(str(species_genome.split('/')[-1]).split('_')[:2])
# Wanted number of threads at the same time:
n_threads = int(sys.argv[3])
# Follow up
follow_up = str(sys.argv[4])


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
# Extract ranges as a list of ranges
###
def extract_ranges(ranges_file_path):
    ranges = list()
    with open(ranges_file_path, 'r') as ranges_file:
        for each_line in ranges_file:
            ranges.append([int(each_line.strip().split('\t')[0]), int(each_line.strip().split('\t')[1])])
    return ranges


###
# Using all the precomputed ranges of factors, will find the pieces in between
###
def extract_uncategorized(record, all_files, species_proxies_directory):
    # This list will contain all the record ranges
    record_ranges = list()

    record_files = [each_file for each_file in all_files if os.path.basename(each_file) == record.id]
    for each_file in record_files:
        factor_ranges = extract_ranges(each_file)
        record_ranges.extend(factor_ranges)
    record_ranges.sort()

    # Merge two intervals which are next to each other
    length_of_list = len(record_ranges)

    i = 0
    while i < (length_of_list - 1):
        if record_ranges[i][1] in (record_ranges[i + 1][0], record_ranges[i + 1][0] - 1):
            record_ranges[i:i + 2] = [[record_ranges[i][0], record_ranges[i + 1][1]]]
            length_of_list -= 1
        else:
            i += 1

    very_first_range = record_ranges[0]

    # We will remove each time the first range to compare to the next, until we only have one range left
    uncategorized = list()
    while len(record_ranges) > 1:
        first_range = record_ranges.pop(0)
        uncategorized.append([first_range[1] + 1, record_ranges[0][0] - 1])

    # if the first start is not the start of the record, must add one first range
    if very_first_range[0] != 0:
        uncategorized.insert(0, [0, very_first_range[0] - 1])

    # On the other end, if the last range is not until the end, we need to add this last range
    if record_ranges[0][1] != len(record):
        uncategorized.append([record_ranges[0][1] + 1, len(record)])



    # If we find some ranges -> can write them
    if uncategorized:
        uncategorized_directory = '/'.join([species_proxies_directory, 'uncategorized', record.id])
        checking_parent(uncategorized_directory)

        # We will now store the clean ranges
        with open(uncategorized_directory, 'w') as outfile:
            for each_range in uncategorized:
                outfile.write('\t'.join([str(each) for each in each_range]) + '\n')


# Fetch this species records for the lengths
records = fetch_fasta(species_genome)

# Directory which lead to the different species proxies
species_proxies_directory = '/'.join(["../files/factor_proxies", species])

# Find all the different factors ready
all_factors = [each for each in extract_path(species_proxies_directory + '/', '*')]
all_files = list(itertools.chain(*[extract_path(each + '/', '*') for each in all_factors]))

Parallel(n_jobs=n_threads)(delayed(extract_uncategorized)
                           (records[each_record], all_files, species_proxies_directory)
                           for each_record in range(len(records)))

# Follow the progression of the analysis
checking_parent(follow_up)
with open(follow_up, 'w') as file:
    file.write('')
