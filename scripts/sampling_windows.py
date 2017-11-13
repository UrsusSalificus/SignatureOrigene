#!/usr/bin/env python3

""" Extract samples windows out of the species genome
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import math
import numpy as np

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Species genome path:
species_genome = str(sys.argv[1])
# Number of samples
n_samples = int(sys.argv[2])
# Wanted window size:
window_size = int(sys.argv[3])
# Output path:
outfile = str(sys.argv[4])


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
# Sample among the genome the wanted number of windows
# Inputs:
#   - records : fetched sequence (fasta)
#   - n_samples : number of sample windows
#   - window_size : size of the wanted sample windows
#   - outfile : path ot the output file
# Output:
#   - A fasta file containing all the windows in separate records
###
def sample_windows(records, n_samples, window_size, outfile):
    # How many windows can we fit into this genome
    max_number_windows = int(sum([math.floor(len(each.seq) / window_size) for each in records]))

    # This list will be fill by the window sample records
    all_sample_records = list()
    # Check if big enough to have the wanted number of sample windows
    if max_number_windows > n_samples:
        # Randomly sampling among all the possible windows of the whole genome
        sample_windows = np.random.randint(0, max_number_windows, size=n_samples)
        sample_windows.sort()
        i = 0  # keep track of at which window sample we are
        sum_n_windows = 0  # keep track of the number of windows we went through
        # Find the window which were randomly sampled in each record
        for each_record in range(len(records)):
            record_length = len(records[each_record].seq)
            record_id = records[each_record].id
            record_n_windows = math.floor(record_length / window_size)
            # While the sample windows are in this record
            while sample_windows[i] < record_n_windows + sum_n_windows:
                start = (sample_windows[i] - sum_n_windows) * window_size
                window_sample_seq = records[each_record].seq[start:start + window_size]
                window_sample_record = SeqRecord(seq=window_sample_seq, id=record_id + '_' + str(start))
                all_sample_records.append(window_sample_record)
                i += 1
                if i == n_samples:
                    break
            sum_n_windows += record_n_windows
    # Else, we will take all the available windows
    else:
        for each_record in range(len(records)):
            n_windows = int(math.floor(len(records[each_record].seq) / window_size))
            names = [records[each_record].id] * n_windows
            windows = np.arange(0, n_windows * window_size, window_size)
            ids = [m + '_' + str(n) for m, n in zip(names, windows)]
            for each_sample in ids:
                start = int(each_sample.split('_')[2])
                window_sample_seq = records[each_record].seq[start:start + window_size]
                window_sample_record = SeqRecord(seq=window_sample_seq, id=each_sample)
                all_sample_records.append(window_sample_record)

    # Write these new records of window sample as a fasta file
    SeqIO.write(all_sample_records, outfile, "fasta")


# Fetch the whole genome records
records = fetch_fasta(species_genome)

sample_windows(records, n_samples, window_size, outfile)
