#!/usr/bin/env python3

"""Extract the percentage of windows' nucleotides which are Low Complexity Regions (LCR)
"""
from Bio import SeqIO
import math
import sys
import os

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Species genome path:
species_genome = str(sys.argv[1])
# Species abbreviation:
species = '_'.join(str(species_genome.split('/')[-1]).split('_')[:2])
# Species feature table path:
species_table = str(sys.argv[2])
# Wanted window size:
window_size = int(sys.argv[3])
# Wanted feature:
feature = str(sys.argv[4])
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


def extract_features (records, species_table, output, id_column, start_column, end_column):
    # Checking parent directory of output are present
    checking_parent(output)

    with open(species_table, 'r') as feature_table, open(output, 'w') as outfile:
        # Must skip the header, which is different among different feature:
        feature_table.readline()
        feature_table.readline()
        feature_table.readline()
        # We will keep this start point:
        start_table = feature_table.tell()

        # From now on, will move through the document by .readline()
        # Problem, tables have different separation type:
        actual_line = feature_table.readline().rsplit()

        for each_record in range(len(records)):
            # Each element of this list represents a nucleotide
            record_proxy = [0] * len(records[each_record].seq)

            # Whenever we are not already at our chromosome part -> skip until at it
            while records[each_record].id != actual_line[id_column]:
                actual_line = feature_table.readline().split()
                # This trick works for most genome, but for A. thaliana genome, we had to take out some
                # records which were at the beginning (of the records order), which makes it scan the whole file
                # In this case, we have to rescan everything (TODO: find a better way to do this)
                try:
                    actual_line[1]
                except:
                    feature_table.seek(start_table)
                    actual_line = feature_table.readline().rsplit()

            # While repeats from the actual record, continue extracting
            try:
                while records[each_record].id == actual_line[id_column]:
                    # If it is the wanted feature
                    if actual_line[feature_column].split('/')[0].strip('?') in feature_type:
                        # For each nucleotide from start to end of the repeat:
                        for each_nucleotide in range(int(actual_line[start_column]), int(actual_line[end_column])):
                            record_proxy[each_nucleotide] = 1
                    actual_line = feature_table.readline().rsplit()
                    # If we get at the last line, actual_line only have one empty entry, which can be detected by
                    # calling the second element ([1])
                    try:
                        actual_line[1]
                    except:
                        break
                if len(records[each_record].seq) > window_size:
                    n_windows = math.floor(len(records[each_record].seq) / window_size)  # Number of windows
                    for start in range(n_windows):
                        window = str(records[each_record].seq)[(start * window_size):((start * window_size) + window_size)]
                        # If any character in the sequence is NOT a standard nucleotides (including unknown nucleotides)
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

# Particularities of the RepeatMasker output (e.g. which column (the numbers) contains what information):
id_column = 4
feature_column = 10
# The feature type depends on the wanted feature
if feature == 'LCR':
    feature_type = 'Low_complexity'
elif feature == 'TE':
    feature_type = ['DNA', 'LINE', 'LTR', 'SINE', 'Retroposon']
elif feature == 'tandem':
    feature_type = ['Satellite', 'Simple_repeat']
start_column = 5
end_column = 6

extract_features(records, species_table, output, id_column, start_column, end_column)

