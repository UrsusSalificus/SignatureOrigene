#!/usr/bin/env python3

"""Compute the recombination rate of each window for Drosophila Melanogaster
"""
from Bio import SeqIO
import math
import sys
import subprocess

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Species genome path:
species_genome = str(sys.argv[1])
# Wanted window size:
window_size = int(sys.argv[2])


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


# RRC, the tool used to calculate the recombination rate of the window, has specific boundaries of the
# chromosomes from the 5.36 version of FlyBase.
# Fetch the 5.36 version
records = fetch_fasta(species_genome)

# Recombination rate is available only on certain records
available = ['2L', '2R', '3L', '3R', 'X']

all_windows = '../tools/RRC/all_windows'
with open(all_windows, 'w') as outfile:
    for each_record in range(len(records)):
        if records[each_record].id in available and len(records[each_record].seq) > window_size:
            n_windows = math.floor(len(records[each_record].seq) / window_size)  # Number of windows
            for each_window in range(n_windows):
                start = each_window * window_size
                end = start + window_size
                # If any character in the sequence is NOT a standard nucleotides (including unknown nucleotides),
                # do NOT compute:
                if not any([c not in 'ATCGatcg' for c in records[each_record].seq[start:end]]):
                    to_write = records[each_record].id + ':' + str(start) + '..' + str(end)
                    outfile.write(to_write + '\n')


subprocess.call(['bash', '../scripts/RRC.sh'])




