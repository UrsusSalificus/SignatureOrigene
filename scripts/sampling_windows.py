#!/usr/bin/env python3

""" Extract samples windows out of the species genome
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
import math
import numpy as np
import re

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Species genome path:
species_genome = str(sys.argv[1])
# Number of samples
n_samples = int(sys.argv[2])
# Wanted window size:
window_size = int(sys.argv[3])
# Type of analysis:
analysis = str(sys.argv[4])
# Output path:
output = str(sys.argv[5])


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
# Given a list of start of windows, will fetch the sequences
# Inputs:
#   - records : fetched sequence (fasta)
#   - window_size : size of the wanted sample windows
#   - sample_windows : Numpy array of start of windows randomly sampled
#   - analysis : type fo analysis currently doing
# Output:
#   - A list of Bio.SeqRecord, as long as the specific window only contains ACTG
#   - May thus yield an empty list if all the sample_windows point to windows with N for example
###
def find_right_sample(records, window_size, sample_windows, analysis):
    # This list will be fill by the window sample records
    all_sample_records = list()

    i = 0  # keep track of at which window sample we are
    sum_n_windows = 0  # keep track of the number of windows we went through
    # Find the window which were randomly sampled in each record
    for each_record in range(len(records)):
        record_length = len(records[each_record].seq)
        record_id = records[each_record].id
        record_n_windows = math.floor(record_length / window_size)
        # While the sample windows are in this record
        # We will catch any "index out of range" error -> end of document = break the while
        try:
            while sample_windows[i] < record_n_windows + sum_n_windows:
                start = (sample_windows[i] - sum_n_windows) * window_size
                # If scaling, need both the record name + the start of the window for the factor extraction
                if analysis == 'scaling':
                    sample_id = record_id + '_' + str(start)
                # Otherwise, the record id will be the name of the factor anyway -> no difference anymore
                else:
                    sample_id = record_id
                window_sample_seq = records[each_record].seq[start:start + window_size]
                # We must make sure there is only ATCG in this sequence:
                if not any([c not in 'ATCGatcg' for c in window_sample_seq]):
                    window_sample_record = SeqRecord(seq=window_sample_seq, id=sample_id)
                    all_sample_records.append(window_sample_record)
                i += 1
        except IndexError:
            break
        sum_n_windows += record_n_windows
    return all_sample_records


###
# Sample among the genome the wanted number of windows
# Inputs:
#   - records : fetched sequence (fasta)
#   - n_samples : number of sample windows
#   - window_size : size of the wanted sample windows
# Output:
#   - A SeqRecord full of the sample windows
###
def sample_windows(records, n_samples, window_size):
    # How many windows can we fit into this genome
    max_number_windows = int(sum([math.floor(len(each.seq) / window_size) for each in records]))

    # Check if big enough to have the wanted number of sample windows
    if max_number_windows > n_samples:
        # Set of all available windows
        all_windows_set = np.array(range(max_number_windows))

        # Randomly sampling among all the possible windows of the whole genome
        sample_windows = np.random.choice(all_windows_set, size=n_samples, replace=False)
        sample_windows.sort()

        # We will remove these sampled windows from the available windows
        # Easy case, sample_windows are directly the indexes
        all_windows_set = np.delete(all_windows_set, sample_windows)

        all_sample_records = find_right_sample(records, window_size, sample_windows, analysis)
        to_resample = n_samples - len(all_sample_records)

        while to_resample != 0:
            # As we delete windows with unknown nucleotides, we may end up with an empty set of remaining windows
            if len(all_windows_set) >= to_resample:
                resample_windows = np.random.choice(all_windows_set, size=to_resample, replace=False)
                resample_windows.sort()
            # If it is the case, we must keep what we already have and return the incomplete records
            else:
                break

            # We will remove these re-sampled windows from the available windows
            # More complicated case: must now find to which index the resample_windows correspond to
            for each in resample_windows:
                all_windows_set = np.delete(all_windows_set, np.where(all_windows_set == each))

            resample_records = find_right_sample(records, window_size, resample_windows, analysis)
            if resample_records:
                for each_resample in resample_records:
                    all_sample_records.append(each_resample)
                    to_resample = n_samples - len(all_sample_records)

        # If we are scaling, we will have to match factor later on to each windows = must sort them
        if analysis == 'scaling':
            # Must sort the records to be in the right order of both id and starting position
            all_sample_ids = set()
            for each_sample in all_sample_records:
                sample_id = '_'.join([each_sample.id.split('_')[0], each_sample.id.split('_')[1]])
                all_sample_ids.add(sample_id)
            right_order_ids = [records[each].id for each in range(len(records)) if records[each].id in all_sample_ids]

            ordered_all_samples_records = list()
            for each_id in right_order_ids:
                id_records = list()
                for each_sample in all_sample_records:
                    sample_id = '_'.join([each_sample.id.split('_')[0], each_sample.id.split('_')[1]])
                    if sample_id == each_id:
                        id_records.append(each_sample)
                id_records.sort(key=lambda r: int(r.id.split('_')[2]))
                for each_record in range(len(id_records)):
                    ordered_all_samples_records.append(id_records[each_record])

            all_sample_records = ordered_all_samples_records

    # Else, we will take all the available windows
    else:
        all_sample_records = list()
        for each_record in range(len(records)):
            # How many window for this specific record?
            n_windows = int(math.floor(len(records[each_record].seq) / window_size))
            # Now we need id fo the record to be: record_id_start_of_the_window
            names = [records[each_record].id] * n_windows
            windows = np.arange(0, n_windows * window_size, window_size)
            # Paste the two
            ids = [m + '_' + str(n) for m, n in zip(names, windows)]
            for each_sample in ids:
                # The last element is the start of the window
                start = int(each_sample.split('_')[-1])
                window_sample_seq = records[each_record].seq[start:start + window_size]
                # If using scaling, must be sure to keep the record as it is + not unknown nucleotides
                if analysis == 'scaling':
                    if not any([c not in 'ATCGatcg' for c in window_sample_seq]):
                        window_sample_record = SeqRecord(seq=window_sample_seq, id=each_sample)
                        all_sample_records.append(window_sample_record)
                # Else we already created many random new k-mer by copy-pasting:
                # creating some more by removing N is not a problem anymore
                else:
                    window_sample_seq = re.sub('[^ATCGatcg]', '', str(window_sample_seq))
                    window_sample_seq = Seq(window_sample_seq)
                    sample_id = each_sample.split('_')[0]
                    window_sample_record = SeqRecord(seq=window_sample_seq, id=sample_id)
                    all_sample_records.append(window_sample_record)

    return all_sample_records


# Fetch the whole genome records
records = fetch_fasta(species_genome)

all_sample_records = sample_windows(records, n_samples, window_size)

# Write these new records of window sample as a fasta file
SeqIO.write(all_sample_records, output, "fasta")