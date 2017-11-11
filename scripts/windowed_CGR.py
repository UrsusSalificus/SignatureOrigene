#!/usr/bin/env python3

"""This script compute CGR of all windows of a certain size, in all species
"""
from Bio import SeqIO
import math
import sys
import os
from joblib import Parallel, delayed

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Species genome path:
species_genome = str(sys.argv[1])
# Species abbreviation:
species = '_'.join(str(species_genome.split('/')[-1]).split('_')[:2])
# Wanted window size:
window_size = int(sys.argv[2])
# Wanted number of threads at the same time:
n_threads = int(sys.argv[3])
# Tracking file:
follow_up = str(sys.argv[4])

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
# Compute the Chaos Game Representation (CGR) of a sequence
# Inputs:
#   - records : fetched sequence (fasta) one wants the CGR computed on
#   - outfile : path to the output file, which will contain the x/y CGR coordinates
#           Note: if empty, will return the coordinates instead of writing a file.
# Output:
#   - Either a file, where each line contain a set of x/y coordinates (separated by \t)
#   - Or the coordinates stocked as [[x coordinates], [y coordinates]]
# Performance :
#   - Takes around 1 second for a 300'000 base pairs long sequence.
###
def CGR_coordinates(records, outfile):
    # Prepare the list of x and y coordinates, already at the right size as we will us numeric range in our loop
    xcord = [0] * len(records)
    ycord = [0] * len(records)
    # actual variables are the pointer variables (were are we?), it will change at each nucleotide.
    actual_x = 0.5
    actual_y = 0.5
    # For each nucleotide, from your actual position in the picture, go half-way to this nucleotide corner
    # Nucleotide corner (x axis, y axis): A (0,0), T (1,0), C (0,1) and G (1,1)
    for each_char in range(len(records)):
        if records[each_char] == "A" or records[each_char] == "a":
            xcord[each_char] = actual_x + ((-actual_x) / 2)
            ycord[each_char] = actual_y + ((-actual_y) / 2)
        elif records[each_char] == "C" or records[each_char] == "c":
            xcord[each_char] = actual_x + ((-actual_x) / 2)
            ycord[each_char] = actual_y + ((1 - actual_y) / 2)
        elif records[each_char] == "G" or records[each_char] == "g":
            xcord[each_char] = actual_x + ((1 - actual_x) / 2)
            ycord[each_char] = actual_y + ((1 - actual_y) / 2)
        elif records[each_char] == "T" or records[each_char] == "t":
            xcord[each_char] = actual_x + ((1 - actual_x) / 2)
            ycord[each_char] = actual_y + ((-actual_y) / 2)
        else:
            raise ValueError('Nucleotide not valid (not ATCG), may impede the rest of the analysis, please clean the '
                             'sequence')
        actual_x = xcord[each_char]
        actual_y = ycord[each_char]

    # If outfile is non-empty, write the output
    if outfile:
        checking_parent(outfile)
        with open(outfile, 'w') as file:
            for each_char in range(len(records)):
                file.write(str(xcord[each_char]) + '\t')
                file.write(str(ycord[each_char]) + '\n')

    # If no 2nd argument was given, outfile is empty (= considered False)
    else:
        # Store the CGR in a form of a list of list
        coordinates = [xcord, ycord]
        return coordinates


###
# Slight update to make the computing CGR function to ignore sequences containing unknown nucleotides
###
def N_sensitive_CGR(window, outfile):
    # If any character in the sequence is NOT a standard nucleotides (including unknown nucleotides), do NOT compute:
    if not any([c not in 'ATCGatcg' for c in window]):
        CGR_coordinates(window, outfile)

# Fetching the genomic fasta file
records = fetch_fasta(species_genome)

# We will know compute the CGR of all windows, in all records:
for each_record in range(len(records)):
    if len(records[each_record].seq) > window_size:
        n_windows = math.floor(len(records[each_record].seq) / window_size)  # Number of windows

        # The follow_up file enable us to know if we are working on "scaling" or "masking/purifying",
        # and will change where we store the CGRs:
        if follow_up.split('/')[0] == '..':
            # We are in the scaling case, where we use the genome "as it is"
            # We store CGRs in source directory:
            seq_directory = '/'.join(['../files/CGRs', str(window_size), species, records[each_record].id, 'CGR_region_'])
        else:
            # We are in the masking/purifying case, where we used sequence of the factor only
            # We store CGRs directly in the purifying directory
            factor = follow_up.split('_')[2]
            seq_directory = '/'.join(['files/CGRs', str(window_size), species, factor,
                                      records[each_record].id, 'CGR_region_'])

        # Parallel the CGR on n_jobs core
        Parallel(n_jobs=n_threads)(delayed(N_sensitive_CGR)
                                   (str(records[each_record].seq[
                                        (start * window_size):((start * window_size) + window_size)]),
                                    seq_directory + str(start))
                                   for start in range(0, n_windows))

# Follow the progression of the analysis
checking_parent(follow_up)
with open(follow_up, 'w') as file:
    file.write('')
