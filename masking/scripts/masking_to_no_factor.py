#!/usr/bin/env python3

""" Masking to have nucleotides NOT composed of the wanted factor
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
import os
import numpy as np
import math
import itertools
import os.path

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Wanted factor:
factor = str(sys.argv[1])
# Wanted window size:
window_size = int(sys.argv[2])
# Sample size:
n_samples = int(sys.argv[3])
# Species genome path:
species_genome = str(sys.argv[4])
# Species abbreviation:
species = '_'.join(str(species_genome.split('/')[-1]).split('_')[:2])
# Output file path
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
# Translate a list of integers into a list of all the ranges found in this list of integers
###
def as_ranges(list_of_integers):
    for p in itertools.groupby(enumerate(list_of_integers), lambda x_y: x_y[1]-x_y[0]):
        b = list(p[1])
        yield b[0][1], b[-1][1]


###
# Build a numpy proxy record of all the indexes where there is NOT the wanted factor, using the record ranges
###
def build_proxy(record_ranges, masked_only_length):
    # This list will contain all the indexes of nucleotide that are within a factor
    record_proxy = np.zeros(masked_only_length, dtype=np.int64)

    end_previous_range = 0
    keep_track = 0
    with open(record_ranges, 'r') as range_file:
        for each_range in range_file:
            line_range = [int(each_range.split()[0]), int(each_range.split()[1]) + 1]

            if end_previous_range == 0:
                start = 0
            else:
                start = end_previous_range
            end = line_range[0] - 1
            non_factor_range = range(start, end)

            # For each nucleotide from start to end of the factor:
            for index, item in enumerate(non_factor_range):
                record_proxy[index + keep_track] = item
            end_previous_range = line_range[1]
            keep_track += len(non_factor_range)

    # We probably don't have factors ranging up to the end of the record sequence -> must add this last range
    if keep_track < masked_only_length:
        last_range = range(keep_track, masked_only_length)
        for index, item in enumerate(last_range):
            record_proxy[index + keep_track] = item

    return record_proxy


###
# Count all the masked factor record length
###
def count_masked_record_length(record, proxies_directory):
    record_ranges = proxies_directory + '/' + record.id

    length_masked_factor = 0
    try:
        with open(record_ranges, 'r') as all_ranges:
            for each_range in all_ranges:
                line_range = each_range.strip().split()
                # Add the range length
                length_masked_factor += len(range(int(line_range[0]), int(line_range[1]) + 1))

            # We want the number of nucleotides which are NOT within the factor
            length_masked_factor = len(record.seq) - length_masked_factor

            return length_masked_factor
    # If no record for this factor, return length of 0
    except FileNotFoundError:
        return 0


###
# Find all the coordinates ranges of this sample of NOT factor
###
def find_sample_ranges_mask(record_ranges, start, window_size):
    sample_ranges = list()

    # This vector will first help us find which ranges are the sample ranges
    adding_up = 0
    with open(record_ranges, 'r') as range_file:
        range = range_file.readline().strip().split()
        # This will be the proxy of start for now, but later help find NOT factor
        adding_up += int(range[0])
        while adding_up < start:
            previous_range = range
            range = range_file.readline().strip().split()
            adding_up += (int(range[0]) - 1) - (int(previous_range[1]) + 1)

        # If we went too far, we must cut the actual range a bit, to find the starting indexes
        if adding_up > start:
            # We use another filler to find when we have a complete window
            left = adding_up - start
            sample_ranges.append([(int(range[0]) - 1) - left, (int(range[0]) - 1)])
        # Otherwise we only have the first nucleotide
        else:
            left = 1
            sample_ranges.append([int(range[1]), int(range[1])])

        while left < window_size:
            previous_range = range
            range = range_file.readline().strip().split()
            left += (int(range[0]) - 1) - (int(previous_range[1]) + 1)
            sample_ranges.append([(int(previous_range[1]) + 1), (int(range[0]) - 1)])

        # Same for the last one, except if left = window size -> nothing to do, perfect
        if left > window_size:
            last_range = sample_ranges.pop()
            right = left - window_size
            sample_ranges.append([last_range[0], last_range[1] - right])

    return sample_ranges


###
# Given a list of start of windows, will fetch the sequences
# Inputs:
#   - records : fetched sequence (fasta)
#   - window_size : size of the wanted sample windows
#   - sample_windows : Numpy array of start of windows randomly sampled
#   - factor_record_lengths : length of all the masked factor records
# Output:
#   - A list of Bio.SeqRecord, as long as the specific window only contains ACTG
#   - May thus yield an empty list if all the sample_windows point to windows with N for example
###
def find_right_sample(records, window_size, sample_windows, factor_record_lengths):
    # This list will be fill by the window sample records
    all_sample_records = list()

    i = 0  # keep track of at which window sample we are
    sum_n_windows = 0  # keep track of the number of windows we went through
    # Find the window which were randomly sampled in each record
    for each_record in range(len(records)):
        record = records[each_record]
        record_ranges = proxies_directory + '/' + record.id

        # Does the rest only if the record exists
        if os.path.isfile(record_ranges):
            # Extract the indexes of nucleotide within the factor
            masked_only_length = factor_record_lengths[each_record]

            record_n_windows = math.floor(masked_only_length / window_size)
            # We will catch any "index out of range" error -> end of document = break the while
            try:
                # While the sample windows are in this record
                while sample_windows[i] < record_n_windows + sum_n_windows:
                    start = (sample_windows[i] - sum_n_windows) * window_size

                    # Find the right index ranges of this sample window
                    sample_factor_only_ranges = find_sample_ranges_mask(record_ranges, start, window_size)

                    # For each of these index ranges -> find the nucleotide associated with
                    sample_seq = str()
                    for each_range in sample_factor_only_ranges:
                        sample_seq += record.seq[each_range[0]:each_range[1] + 1]

                    # We must make sure there is only ATCG in this sequence, else simply not take this sample
                    if not any([c not in 'ATCGatcg' for c in sample_seq]):
                        window_sample_record = SeqRecord(seq=sample_seq, id=factor)
                        all_sample_records.append(window_sample_record)

                    # Moving to the next sample
                    i += 1
            except IndexError:
                break
            sum_n_windows += record_n_windows
    return all_sample_records


###
# Extract from the proxy all the nucleotide that are from the wanted factor
# Inputs:
#   - records : fetched sequence (fasta)
#   - proxies_directory : directory of the proxies obtained through the factor_proxies script
#   - window_size : size of the wanted sample windows
#   - n_samples : wantted number of sample windows
# Output:
#   - A list of n samples as SeqRecords of masked factor
###
def sampling_using_proxies(records, proxies_directory, window_size, n_samples):
    # Find all the masked record total lengths, before any sampling
    factor_record_lengths = list()
    for each_record in range(len(records)):
        factor_record_lengths.append(count_masked_record_length(records[each_record], proxies_directory))

    # Find the maximum number of sample windows we can get out of these masked record
    max_number_windows = int(sum([math.floor(each / window_size) for each in factor_record_lengths]))

    # Check if big enough to have the wanted number of sample windows
    # But if really big = some memory issues with these big sequences of masked factor
    # As such, we will not create it, but use the ranges to sample inside this sequences of masked factor
    if max_number_windows > n_samples:
        # Set of all available windows
        all_windows_set = np.array(range(max_number_windows))

        # Randomly sampling among all the possible windows of the whole genome
        sample_windows = np.random.choice(all_windows_set, size=n_samples, replace=False)
        sample_windows.sort()

        # We will remove these sampled windows from the available windows
        # Easy case, sample_windows are directly the indexes
        all_windows_set = np.delete(all_windows_set, sample_windows)

        # Use the find_right_sample to extract the right sequence using this set of sample windows
        all_sample_records = find_right_sample(records, window_size, sample_windows, factor_record_lengths)
        # After finding the records, we must check if we did not lost some due to unknown nucleotide
        to_resample = n_samples - len(all_sample_records)

        # If it is the case, we must try while we still have some sample windows at hand
        while to_resample != 0:
            # As we delete windows with unknown nucleotides, we may end up with an empty set of remaining windows
            if len(all_windows_set) >= to_resample:
                resample_windows = np.random.choice(all_windows_set, size=to_resample, replace=False)
                resample_windows.sort()
            # If it is the case, we must keep what we already have and return the incomplete sample records
            else:
                break

            # We will remove these re-sampled windows from the available windows
            # More complicated case: must now find to which index the resample_windows correspond to
            for each in resample_windows:
                all_windows_set = np.delete(all_windows_set, np.where(all_windows_set == each))

            resample_records = find_right_sample(records, window_size, resample_windows, factor_record_lengths)
            # If did found some new sample windows, append it to the main list, else continue the search
            if resample_records:
                for each_resample in resample_records:
                    all_sample_records.append(each_resample)
                    to_resample = n_samples - len(all_sample_records)
    # Else, the records are quite small, we will thus take all the available windows
    # As size is not a problem anymore, we can easily build the sequences of masked factor
    else:
        # This list will contain all the sample records
        all_sample_records = list()

        # If too small -> we must be able to concatenate everything
        too_small = str()

        # We will process the factor from all records
        for each_record in range(len(records)):
            record = records[each_record]
            record_ranges = proxies_directory + '/' + record.id

            # Does the rest only if the record exists
            if os.path.isfile(record_ranges):
                # We need the indexes of all nucleotide that are within the wanted factor
                masked_only_length = factor_record_lengths[each_record]
                masked_only = build_proxy(record_ranges, masked_only_length)

                # We translate the whole genome sequence of this record as a numpy array to use with the indexes
                numpy_seq = np.asarray([c for c in record.seq])

                # Extract the sequence of nucleotides which are within the wanted factor
                masked_seq = numpy_seq[masked_only]

                # We already created many random new k-mer by copy-pasting:
                # creating some more by removing N is not a problem anymore
                masked_seq = str(''.join([c for c in masked_seq if c in 'ATCGatcg']))

                # Check if it is not too small already:
                if masked_only_length < window_size:
                    too_small += masked_seq
                else:
                    # How many window for this specific record?
                    record_n_windows = int(math.floor(len(masked_seq) / window_size))

                    # For each sample, extract the right nucleotides and append it
                    for each_sample in range(record_n_windows):
                        start = each_sample * window_size

                        sample_seq = Seq(masked_seq[start:start + window_size])

                        all_sample_records.append(SeqRecord(seq=sample_seq, id=factor))

        # Check if we can get something out of too_small:
        if too_small and len(too_small) > window_size:
            # How many window for this concatenated sequence?
            concatenated_n_windows = int(math.floor(len(too_small) / window_size))
            for each_sample in range(concatenated_n_windows):
                start = each_sample * window_size

                concatenated_seq = Seq(too_small[start:start + window_size])

                all_sample_records.append(SeqRecord(seq=concatenated_seq, id=factor))

    return all_sample_records


# Fetch all the records from this species fasta
records = fetch_fasta(species_genome)

# Directory containing all the ranges in all the different files
proxies_directory = '/'.join(['../files/factor_proxies', species, factor])

all_windows_samples = sampling_using_proxies(records, proxies_directory, window_size, n_samples)

# Writing the record as fasta file
SeqIO.write(all_windows_samples, output, "fasta")