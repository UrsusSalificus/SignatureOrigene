#!/usr/bin/env python3

"""Extract the percentage of windows' average recombination rate
"""
from Bio import SeqIO
import math
import sys
import os
import subprocess

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Wanted window size:
window_size = int(sys.argv[1])
# Species sample windows path:
species_sample = str(sys.argv[2])
# Species abbreviation:
species = '_'.join(str(species_sample.split('/')[-1]).split('_')[:2])
# Temporary file containing the sample windows coordinates
species_temp = species + '_temp_all_windows'
# Species recombination rate as spline functions:
species_RR = str(sys.argv[3])
# Output path:
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
    return records


###
# Write a BED file containing the coordinates (first/last) of each windows not containing any N in them
# Input:
#   - records : fetched sequence (fasta) one wants the CGR computed on
#   - window_size : size of the sliding window
#   - species_temp : path of the file which will be written on (removed at the end of the script)
# Output:
#   - A BED file, each row corresponds to a window, and on each row one can find the record ids and
#       boundaries of the window.
###
def all_coordinates_BED(records, window_size, species_temp):
    with open(species_temp, 'w') as outfile:
        for each_record in range(len(records)):
            window_id = '_'.join([records[each_record].id.split('_')[0], records[each_record].id.split('_')[1]])
            start = int(records[each_record].id.split('_')[2])
            end = start + window_size
            # If any character in the sequence is NOT a standard nucleotides (including unknown nucleotides),
            # do NOT compute:
            if not any([c not in 'ATCGatcg' for c in records[each_record].seq[start:end]]):
                to_write = '\t'.join([window_id, str(start), str(end)])
                outfile.write(to_write + '\n')


# Fetch all the records from this species fasta
samples = fetch_fasta(species_sample)

all_coordinates_BED(samples, window_size, species_temp)

# Actual averaging of the spline fit of the recombination rates, for each windows
subprocess.call(['Rscript', 'scripts/extract_RR.R', species_temp, species_RR, output])

# Will remove the temp file
os.remove(species_temp)
