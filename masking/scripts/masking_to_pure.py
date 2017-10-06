#!/usr/bin/env python3

""" Masking to only have nucleotides composed of the wanted feature
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
# Species genome path:
species_genome = str(sys.argv[2])
# Species feature table path, depends on the type of factor:
if factor in ['LCR', 'TE', 'tandem']:
    species_table = str(sys.argv[3])
elif factor in ['CDS', 'RNA']:
    species_table = str(sys.argv[4])
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


def reading_line(factor, feature_table):
    if factor in ['LCR', 'TE', 'tandem']:
        return feature_table.readline().rsplit()
    elif factor in ['CDS', 'RNA']:
        return feature_table.readline().split('\t')


def True_if_right_factor_strand(factor, actual_line, feature_column, feature_type, strand_column):
    if factor in ['LCR', 'TE', 'tandem']:
        return actual_line[feature_column].split('/')[0].strip('?') in feature_type
    elif factor in ['CDS', 'RNA']:
        return actual_line[feature_column] in feature_type and actual_line[strand_column] == '+'


def extract_factor(records, factor, species_table, output, id_column, feature_column, feature_type, strand_column,
                start_column, end_column):
    # Checking parent directory of output are present
    checking_parent(output)

    # We will store the records in this list
    factor_records = list()

    with open(species_table, 'r') as feature_table, open(output, 'w') as outfile:
        # Must skip the header (which differs in between feature_table and repeats:
        if factor in ['LCR', 'TE', 'tandem']:
            feature_table.readline()
            feature_table.readline()
            feature_table.readline()
        elif factor in ['CDS', 'RNA']:
            feature_table.readline()

        # From now on, will move through the document by .readline()
        actual_line = reading_line(factor, feature_table)

        for each_record in range(len(records)):
            # This string will contain all the nucleotide of the record which are coding
            all_ranges = list()

            # Whenever we are not already at our chromosome part -> skip until at it
            while records[each_record].id != actual_line[id_column]:
                actual_line = reading_line(factor, feature_table)

            # We also have to find the first time the wanted feature appears
            while not True_if_right_factor_strand(factor, actual_line, feature_column, feature_type, strand_column):
                actual_line = reading_line(factor, feature_table)

            # This line will be the first result
            all_ranges.append([int(actual_line[start_column]),int(actual_line[end_column])])

            # It now become the precedent line
            precedent_line = actual_line
            actual_line = reading_line(factor, feature_table)

            # While from the actual record, continue extracting
            while records[each_record].id == actual_line[id_column]:
                # Only do this for wanted feature
                if True_if_right_factor_strand(factor, actual_line, feature_column, feature_type, strand_column):
                    # We must check if there is an overlap
                    # 1) If there is overlap, we must cut the intersection out
                    if int(actual_line[start_column]) <= int(precedent_line[end_column]):
                        # 1.2) Now, it may be that a small CDS may be inside a bigger CDS (thus passing the if).
                        # We will check for this and only store the CDS range if it is not the case
                        if int(actual_line[end_column]) > int(precedent_line[end_column]):
                            precedent_range = set(range(int(precedent_line[start_column]),
                                                        int(precedent_line[end_column])))
                            actual_range = set(range(int(actual_line[start_column]), int(actual_line[end_column])))
                            upper_part = actual_range - precedent_range
                            all_ranges.append([min(upper_part), max(upper_part)])
                    # 2) If no overlap, simply store the CDS range
                    else:
                        all_ranges.append([int(actual_line[start_column]), int(actual_line[end_column])])
                    # It now become the precedent line
                    precedent_line = actual_line
                # Continue to next line
                actual_line = reading_line(factor, feature_table)
                # If we get at the last line, actual_line only have one empty entry, which can be detected by
                # calling the second element ([1])
                try:
                    actual_line[1]
                except:
                    break

            # We can finally store all the sequences as a single sequence
            factor_only = str()
            for each_range in all_ranges:
                factor_only += records[each_record].seq[each_range[0]:each_range[1]]

            new_record = SeqRecord(seq = factor_only, id = records[each_record].id)
            factor_records.append(new_record)

        # Write the new list of records
        SeqIO.write(factor_records, outfile, "fasta")


if factor in ['LCR', 'TE', 'tandem']:
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
elif factor in ['CDS', 'RNA']:
    id_column = 6
    feature_column = 0
    # The feature type depends on the wanted feature
    if factor == 'CDS':
        feature_type = 'CDS'
    elif factor == 'RNA':
        feature_type = ['misc_RNA', 'ncRNA', 'rRNA', 'tRNA']
    strand_column = 9
    start_column = 7
    end_column = 8

# Fetch all the records from this species fasta
records = fetch_fasta(species_genome)

extract_factor(records, factor, species_table, output, id_column, feature_column, feature_type, strand_column,
                start_column, end_column)