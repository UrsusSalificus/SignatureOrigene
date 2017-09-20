#!/usr/bin/env python3

"""Compute the recombination rate of each window for Drosophila Melanogaster
"""
from Bio import SeqIO
import math
import sys
import os
import subprocess

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Species genome path:
species_genome = str(sys.argv[1])
# Species recombination rates path:
RR = str(sys.argv[2])
# Wanted window size:
window_size = int(sys.argv[3])
# Output path:
output = str(sys.argv[4])


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


records = fetch_fasta(species_genome)
n_windows = math.floor(len(records[0].seq) / window_size)  # Number of windows

species_genome='../data/genomes/hsap_36sample_genomes.fna'
window_size= 15000

with open(output, 'w') as outfile, open(RR, 'r') as rates:
    # Skip header
    rates.readline()
    actual_line = rates.readline().split()
    position = float(actual_line[0])
    recombination_rate = float(actual_line[1])
    for each_window in range(n_windows):
        all_RR = []
        start = each_window * window_size
        end = start + window_size
        # If any character in the sequence is NOT a standard nucleotides (including unknown nucleotides),
        # do NOT compute:
        if not any([c not in 'ATCGatcg' for c in records[0].seq[start:end]]):
            while position <= end:
                all_RR.append(recombination_rate)
                actual_line = rates.readline().split()
                # If we get at the last line, actual_line only have one empty entry, which can be detected by
                # calling the second element ([1])
                try:
                    actual_line[1]
                except:
                    break
                position = float(actual_line[0])
                recombination_rate = float(actual_line[1])
            if len(all_RR) == 0:
                outfile.write('NA' + '\n')
            else:
                RR_average = sum(all_RR) / len(all_RR)
                outfile.write(str(RR_average) + '\n')







