#!/usr/bin/env python3

""" Extracting percentage of nucleotides that are within the factor
"""
from Bio import SeqIO
import sys
import os
import numpy as np

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Wanted factor:
factor = str(sys.argv[1])
# Wanted window size:
window_size = int(sys.argv[2])
# Species whole genome path:
species_genome = str(sys.argv[3])
# Species abbreviation:
species = '_'.join(str(species_genome.split('/')[-1]).split('_')[:2])
# Species sample windows path:
species_sample = str(sys.argv[4])
# Output path:
output = str(sys.argv[5])


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
# Build a numpy proxy record of factor using a files with the ranges of nucleotide within this factor + record length
###
def build_proxy(record_ranges, record_length):
    # Each element of this list represents a nucleotide
    record_proxy = np.zeros(record_length, dtype=np.int8)

    with open(record_ranges, 'r') as range_file:
        for each_range in range_file:
            start = int(each_range.split()[0])
            end =  int(each_range.split()[1])

            # For each nucleotide from start to end of the CDS:
            for each_nucleotide in range(start, end):
                record_proxy[each_nucleotide] = 1

    record_proxy = np.asarray(record_proxy)

    return record_proxy


###
# Extract the percentage of the sample windows nucleotides which are linked ot the wanted factor
# Inputs:
#   - record_proxy : proxy obtained through the build_proxy function
#   - record_samples : all the sample sequences of this record
#   - window_size : size of the chosen window
# Output:
#   - List containing for each sample, the percentage of wanted factor
###
def extract_percentage (record_proxy, record_samples, window_size):
    # We will store all the percentages of factor per windows
    record_percentages = list()

    for each_sample in record_samples:
        window_start = int(each_sample.id.split('_')[2])
        window_percentage = sum(record_proxy[window_start:window_start + window_size])
        record_percentages.append((window_percentage / window_size) * 100)

    return record_percentages


# Fetch all the records from this species fasta
records = fetch_fasta(species_genome)
samples = fetch_fasta(species_sample)

# Directory containing all the ranges in all the different files
proxies_directory = '/'.join(['../files/factor_proxies', str(window_size), species, factor])

samples_percentages = list()
for each_record in range(len(records)):
    record_length = len(records[each_record].seq)
    record_id = records[each_record].id
    record_ranges = proxies_directory + '/' + record_id
    # Build a proxy of record where 1 = nucleotide within factor ; 0 = nucleotide out of factor
    record_proxy = build_proxy(record_ranges, record_length)

    # Extracting all the samples from this record
    record_samples = [each for each in samples if '_'.join(each.id.split('_')[0:2]) == record_id]

    # Check if not empty
    if record_samples:
        record_percentages = extract_percentage (record_proxy, record_samples, window_size)
        samples_percentages.extend(record_percentages)

# Checking parent directory of output are present
checking_parent(output)

# Writing the fasta file
with open(output, 'w') as outfile:
    for each_sample in range(len(samples_percentages)):
        outfile.write(samples[each_sample].id + '\t' + str(samples_percentages[each_sample]) + '\n')