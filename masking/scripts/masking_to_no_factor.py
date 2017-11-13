#!/usr/bin/env python3

"""This script compute k-mer frequencies using CGR of all windows of a certain size
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Wanted factor:
factor = str(sys.argv[1])
# Wanted window size:
window_size = int(sys.argv[2])
# Species genome path:
species_genome = str(sys.argv[3])
# Species feature table path, depends on the type of factor:
if factor in ['LCR', 'TE', 'tandem']:
    species_table = str(sys.argv[4])
    factor_type = 'repeats'
elif factor in ['CDS', 'RNA', 'intron', 'UTR']:
    species_table = str(sys.argv[5])
    factor_type = 'features'
# Output path:
output = str(sys.argv[6])


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


def reading_line(factor_type, feature_table):
    if factor_type == 'repeats':
        return feature_table.readline().rsplit()
    else:
        return feature_table.readline().split('\t')


def True_if_right_factor_strand(factor_type, actual_line, feature_column, feature_type, strand_column):
    # Watch out for commentaries (thus length 1)
    if len(actual_line) > 1:
        # If bigger than one -> feature line
        if factor_type == 'repeats':
            return actual_line[feature_column].split('/')[0].strip('?') in feature_type
        else:
            return actual_line[feature_column] in feature_type and actual_line[strand_column] == '+'
    else:
        # Else we know it is a commentary = not the right factor...
        return False


def extract_factor(records, factor, factor_type, species_table, output, id_column, feature_column, feature_type,
                   strand_column, start_column, end_column):
    # Checking parent directory of output are present
    checking_parent(output)

    # We will store the records in this list
    factor_records = list()

    with open(species_table, 'r') as feature_table, open(output, 'w') as outfile:
        # Must skip the header (which differs in between feature_table and repeats:
        if factor_type == 'repeats':
            feature_table.readline()
            feature_table.readline()
            feature_table.readline()
            actual_line = reading_line(factor_type, feature_table)
        else:
            line = feature_table.readline()
            # Must skip the headers (varying length)
            while line.startswith('#'):
                line = feature_table.readline()
            actual_line = line.split('\t')

        for each_record in range(len(records)):
            # Each element of this list represents a nucleotide
            record_proxy = [0] * len(records[each_record].seq)
            # This set will contain all the ranges of our wanted factor
            all_ranges = set()

            # Whenever we are not already at our chromosome part -> skip until at it
            while records[each_record].id != actual_line[id_column]:
                actual_line = reading_line(factor_type, feature_table)

            # We also have to find the first time the wanted feature appears
            while not True_if_right_factor_strand(factor_type, actual_line, feature_column, feature_type,
                                                  strand_column):
                actual_line = reading_line(factor_type, feature_table)

            # This line will be the first result
            all_ranges.add((int(actual_line[start_column]), int(actual_line[end_column])))

            # Continue the search
            actual_line = reading_line(factor_type, feature_table)

            # While from the actual record, continue extracting
            while records[each_record].id == actual_line[id_column]:
                # Only do this for wanted feature
                if True_if_right_factor_strand(factor_type, actual_line, feature_column, feature_type, strand_column):
                    all_ranges.add((int(actual_line[start_column]), int(actual_line[end_column])))
                    # Continue searching
                    actual_line = reading_line(factor_type, feature_table)
                # If it is not our factor, just continue the search
                else:
                    actual_line = reading_line(factor_type, feature_table)
                # If we get at the last line, actual_line only have one empty entry
                if not actual_line:
                    break

            for each_range in all_ranges:
                # For each nucleotide from start to end of the CDS:
                for each_nucleotide in range(each_range[0], each_range[1]):
                    record_proxy[each_nucleotide] = 1

            # We will find all non-overlapping ranges of anything but the factors:
            factor_only = str()
            no_overlap_ranges = list()
            i = 0
            length_record = len(record_proxy)
            # Until the end of record
            while True:
                # We will catch any "index out of range" error -> end of document = break the while
                try:
                    start = i + 1 if i != 0 else i
                    while not record_proxy[i]:
                        i += 1
                        # We must check here if we get to the end of the record -> to have
                        # the range from last factor to end
                        if i == length_record:
                            break
                    end = i - 1
                    no_overlap_ranges.append([start, end])
                    # This will raise an error when at the end of record -> except -> break
                    while record_proxy[i]:
                        i += 1
                except IndexError:
                    break

            for each_range in no_overlap_ranges:
                factor_only += records[each_record].seq[each_range[0]:each_range[1]]

            new_record = SeqRecord(seq = factor_only, id = '_'.join([factor, records[each_record].id]))
            factor_records.append(new_record)

        # We might end up with records which are too small
        # In this case, we must concatenate them:
        if any([len(each_record) < window_size for each_record in factor_records]):
            concatenated_seq = str()
            for each_record in factor_records:
                concatenated_seq += each_record
            concatenated_records = SeqRecord(seq = concatenated_seq.seq, id = factor)

            # Write the new list of records
            SeqIO.write(concatenated_records, outfile, "fasta")
        else:
            # Write the new list of records
            SeqIO.write(factor_records, outfile, "fasta")


if factor_type == 'repeats':
    id_column = 4
    feature_column = 10
    # The feature type depends on the wanted feature
    if factor == 'LCR':
        feature_type = 'Low_complexity'
    elif factor == 'TE':
        feature_type = ['DNA', 'LINE', 'LTR', 'SINE', 'Retroposon']
    elif factor == 'tandem':
        feature_type = ['Satellite', 'Simple_repeat']
    strand_column = False   # NOT USED
    start_column = 5
    end_column = 6
else:
    id_column = 0
    feature_column = 2
    # The feature type depends on the wanted feature
    if factor == 'CDS':
        feature_type = 'CDS'
    elif factor == 'RNA':
        feature_type = ['misc_RNA', 'ncRNA', 'rRNA', 'tRNA']
    elif factor == 'intron':
        feature_type = 'intron'
    elif factor == 'UTR':
        feature_type = ['five_prime_UTR', 'three_prime_UTR']
    strand_column = 6
    start_column = 3
    end_column = 4

# Fetch all the records from this species fasta
records = fetch_fasta(species_genome)

extract_factor(records, factor, factor_type, species_table, output, id_column, feature_column, feature_type,
               strand_column, start_column, end_column)