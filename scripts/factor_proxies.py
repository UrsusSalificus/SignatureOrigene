#!/usr/bin/env python3

""" Create a proxy of each records, where each nucleotide within the factor = 1, the rest = 0
"""
from Bio import SeqIO
import sys
import os
import numpy as np
import itertools

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Wanted factor:
factor = str(sys.argv[1])
# Wanted window size:
window_size = int(sys.argv[2])
# Species genome path:
species_genome = str(sys.argv[3])
# Species abbreviation:
species = '_'.join(str(species_genome.split('/')[-1]).split('_')[:2])
# Species feature table path, depends on the type of factor:
if factor in ['LCR', 'TE', 'tandem']:
    species_table = str(sys.argv[4])
    factor_type = 'repeats'
elif factor in ['CDS', 'RNA', 'intron', 'UTR']:
    species_table = str(sys.argv[5])
    factor_type = 'features'
# Tracking file:
follow_up = str(sys.argv[6])


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
# Reading line varies in between species table (repeats vs features)
###
def reading_line(factor_type, feature_table):
    if factor_type == 'repeats':
        return feature_table.readline().rsplit()
    else:
        return feature_table.readline().split('\t')


###
# Will output True if the line contain the right factor (and on the + strand for feature table)
# Will work differently depending on whether working on repeats or features
###
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


###
# Compute a proxy of each records composed of 0 (nucleotide != factor) and 1 (nucleotide == factor)
# Inputs:
#   - records : fetched sequence (fasta) of the species whole genome
#   - factor_type : indicating if wants factor that are either repeats or features
#   - feature_type : what are the pattern to look for in the factor/gene column
#   - species_table : file containing the wanted factor (either RepeatMasker output or gff file)
#   - *_column : various information related to the internal structure of the species_table file
#   - output : path to the output directory
#           Note: if empty, will return the proxy as numpy array
# Output:
#   - Either a numpy array of n (number of different records in records) ranges, if output is non-empty
#       - Each proxy contain m (number of non-overlapping ranges) ranges
#   - Or n (number of records in the fasta file) files containing the non-overlapping ranges
#       of nucleotide which are factors
###
def extract_factor(records, factor_type, feature_type, species_table, id_column, feature_column,
                   strand_column, start_column, end_column, output):
    # Note: we always add -1 to make it compatible the pythonic start of counting at 0!
    if not output:
        # We will store the proxies of records in this list
        proxies_records = list()

    with open(species_table, 'r') as feature_table:
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
            # This set will contain all the ranges of our wanted factor
            all_ranges = set()

            # We ended up with some issues with UTR, as for unknown reasons, the line are sometimes placed
            # after the end of the record...
            if factor == 'UTR':
                for each_line in feature_table:
                    actual_line = each_line.split('\t')
                    if True_if_right_factor_strand(factor_type, actual_line, feature_column, feature_type,strand_column) \
                            and records[each_record].id == actual_line[id_column]:
                        all_ranges.add((int(actual_line[start_column]) - 1, int(actual_line[end_column]) - 1))
            # Any other factor than than UTR does not have this problem -> slightly faster version:
            else:
                # Whenever we are not already at our chromosome part -> skip until at it
                while records[each_record].id != actual_line[id_column]:
                    actual_line = reading_line(factor_type, feature_table)

                # We also have to find the first time the wanted feature appears
                while not True_if_right_factor_strand(factor_type, actual_line, feature_column, feature_type,
                                                      strand_column):
                    actual_line = reading_line(factor_type, feature_table)

                # This line will be the first result
                all_ranges.add((int(actual_line[start_column]) - 1, int(actual_line[end_column]) - 1))

                # Continue the search
                actual_line = reading_line(factor_type, feature_table)

                # While from the actual record, continue extracting
                while records[each_record].id == actual_line[id_column]:
                    # Only do this for wanted feature
                    if True_if_right_factor_strand(factor_type, actual_line, feature_column, feature_type, strand_column):
                        all_ranges.add((int(actual_line[start_column]) - 1, int(actual_line[end_column]) - 1))
                        # Continue searching
                        actual_line = reading_line(factor_type, feature_table)
                    # If it is not our factor, just continue the search
                    else:
                        actual_line = reading_line(factor_type, feature_table)
                    # If we get at the last line, actual_line only have one empty entry
                    if not actual_line:
                        break
            # We will now extract the factor only ranges
            # Each element of this list represents a nucleotide
            record_proxy = [0] * len(records[each_record].seq)

            for each_range in all_ranges:
                # For each nucleotide from start to end of the CDS:
                for each_nucleotide in range(each_range[0], each_range[1]):
                    record_proxy[each_nucleotide] = 1

            # Translate it to numpy array
            record_proxy = np.asarray(record_proxy)

            factor_only = np.where(record_proxy == 1)[0]

            factor_ranges = list(as_ranges(factor_only))

            # If not output, use the ranges to mark the proxy
            if not output:
                # Add the proxy of factor to the list
                proxies_records.append(factor_ranges)
            # Else, write the ranges on a file named as the record id
            else:
                record_file = output + '/' + records[each_record].id

                checking_parent(record_file)
                with open(record_file, 'w') as outfile:
                    for each_range in factor_ranges:
                        outfile.write(str(each_range[0]) + '\t' + str(each_range[1]) + '\n')
    if not output:
        return proxies_records


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
    strand_column = False  # NOT USED
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

proxies_directory = '/'.join(['../files/factor_proxies', str(window_size), species, factor])

# Compute factor proxy of records
extract_factor(records, factor_type, feature_type, species_table, id_column, feature_column,
               strand_column, start_column, end_column, proxies_directory)

# Follow the progression of the analysis
checking_parent(follow_up)
with open(follow_up, 'w') as file:
    file.write('')