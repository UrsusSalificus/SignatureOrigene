#!/usr/bin/env python3

"""Extract the percentage of windows' nucleotides which part of the selected feature
"""
from Bio import SeqIO
import math
import sys
import os

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Wanted window size:
window_size = int(sys.argv[1])
# Wanted feature:
feature = str(sys.argv[2])
# Output path:
output = str(sys.argv[3])
# Species genome path:
species_genome = str(sys.argv[4])
# Species abbreviation:
species = '_'.join(str(species_genome.split('/')[-1]).split('_')[:2])
# Species feature table path:
if feature == "CDS":
    species_table = str(sys.argv[5])
elif feature == "LCR":
    species_table = str(sys.argv[6])


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


def extract_features (feature, species_table, output, id_column, feature_column, feature_type, strand_column,
                      start_column, end_column):
    # Checking parent directory of output are present
    checking_parent(output)

    with open(species_table, 'r') as feature_table, open(output, 'w') as outfile:
        # Must skip the header, which is different among different feature:
        if feature == 'CDS':
            feature_table.readline()
        if feature == 'STR':
            feature_table.readline()
            feature_table.readline()
            feature_table.readline()
        # From now on, will move through the document by .readline()
        actual_line = feature_table.readline().split()
        for each_record in range(len(records)):
            # Each element of this list represents a nucleotide
            record_proxy = [0] * len(records[each_record].seq)
            # While repeats from the actual record, continue extracting
            try :
                while records[each_record].id == actual_line[id_column]:
                    # Now depends on the feature type:
                    if feature == 'CDS':
                        # If it is the wanted feature and on the positive strand
                        if actual_line[feature_column] == feature_type and actual_line[strand_column] == '+':
                            # For each nucleotide from start to end of the CDS:
                            for each_nucleotide in range(int(actual_line[start_column]), int(actual_line[end_column])):
                                record_proxy[each_nucleotide] = 1
                    elif feature == 'LCR':
                        # For each nucleotide from start to end of the repeat:
                        for each_nucleotide in range(int(actual_line[start_column]), int(actual_line[end_column])):
                            record_proxy[each_nucleotide] = 1
                    actual_line = feature_table.readline().split()
                    # If we get at the last line, actual_line only have one empty entry, which can be detected by
                    # calling the second element ([1])
                    try:
                        actual_line[1]
                    except:
                        break
                if len(records[each_record].seq) > window_size:
                    n_windows = math.floor(len(records[each_record].seq) / window_size)  # Number of windows
                    for start in range(0, n_windows):
                        window = str(records[each_record].seq)[(start * window_size):((start * window_size) + window_size)]
                        # If any character in the sequence is NOT a standard nucleotides (including unknown nucleotides),
                        # do NOT compute:
                        if not any([c not in 'ATCGatcg' for c in window]):
                            region_CDS = sum(record_proxy[(start * window_size):((start * window_size) + window_size)])
                            outfile.write(records[each_record].id + '\t')
                            outfile.write(str((region_CDS/window_size)*100) + '\n')
            except:
                print("Extracting stopped at record " + str(each_record) + "/" + str(len(records)))
                print("If near the end, may come from the fact that RepeatMasker was not done on mitochondrion!")
                break

# Fetch all the records from this species fasta
records = fetch_fasta(species_genome)

# Particularities of each feature (e.g. which column (the numbers) contains what information):
if feature == 'CDS':
    id_column = 6
    feature_column = 0
    feature_type = 'CDS'
    strand_column = 9
    start_column = 7
    end_column = 8
elif feature == 'LCR':
    id_column = 4
    feature_column = 10     # NOT USED
    feature_type = ''       # NOT USED
    strand_column = 8       # NOT USED
    start_column = 5
    end_column = 6

extract_features(feature, species_table, output, id_column, feature_column, feature_type, strand_column,
                      start_column, end_column)

