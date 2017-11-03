#!/usr/bin/env python3

""" Extract the percentage of windows' nucleotides which are composed of the wanted feature
(out of the gff file of the species)
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
# Wanted factor:
factor = str(sys.argv[4])
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


def extract_feature(records, species_table, output, id_column, feature_column, feature_type, strand_column,
                start_column, end_column, window_size):
    # Checking parent directory of output are present
    checking_parent(output)

    with open(species_table, 'r') as feature_table, open(output, 'w') as outfile:
        line = feature_table.readline()
        # Must skip the headers (varying length)
        while line.startswith('#'):
            line = feature_table.readline()
        # From now on, will move through the document by .readline()
        actual_line = line.split('\t')
        for each_record in range(len(records)):
            # Each element of this list represents a nucleotide
            record_proxy = [0] * len(records[each_record].seq)

            # Whenever we are not already at our chromosome part -> skip until at it
            while records[each_record].id != actual_line[id_column]:
                actual_line = feature_table.readline().split('\t')

            # While repeats from the actual record, continue extracting
            try:
                while records[each_record].id == actual_line[id_column]:
                    # If it is the wanted feature and on the positive strand
                    if actual_line[feature_column] in feature_type and actual_line[strand_column] == '+':
                        # Here, will now depend whether it is introns we are extracting or any other factor
                        if feature_type == 'exon':
                            # For introns, we need the next line
                            next_line = feature_table.readline().split('\t')
                            # We need the next line to also be exon, indicating an intron between them
                            # We use a while when there is multiple exons
                            while next_line[feature_column] in feature_type:
                                # Mark the nucleotide as factor in the proxy,
                                # by looking at end of first exon to start of the second one
                                for each_nucleotide in range(int(actual_line[end_column])+1,int(next_line[start_column])-1):
                                    record_proxy[each_nucleotide] = 1
                                # Move to the next line
                                actual_line = next_line
                                next_line = feature_table.readline().split('\t')
                        else:
                            # For all other factors, each line containing the right factor is added to the proxy
                            # For each nucleotide from start to end of the factor:
                            for each_nucleotide in range(int(actual_line[start_column]), int(actual_line[end_column])):
                                record_proxy[each_nucleotide] = 1
                    # Continue the search:
                    actual_line = feature_table.readline().split('\t')
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
                            region_introns = sum(record_proxy[(start * window_size):((start * window_size) + window_size)])
                            outfile.write(records[each_record].id + '\t')
                            outfile.write(str((region_introns / window_size) * 100) + '\n')
            except:
                print("Extracting stopped at record " + str(each_record) + "/" + str(len(records)))
                break


# Fetch all the records from this species fasta
records = fetch_fasta(species_genome)

# Particularities of the feature table (e.g. which column contains what information):
id_column = 0
feature_column = 2
# The feature type depends on the wanted feature
if factor == 'CDS':
    feature_type = 'CDS'
elif factor == 'RNA':
    feature_type = ['misc_RNA', 'ncRNA', 'rRNA', 'tRNA']
elif factor == 'intron':
    feature_type = 'exon'
    # Indeed, introns are not directly annotated, but can be infered through exons
elif factor == 'UTR':
    feature_type = ['five_prime_UTR', 'three_prime_UTR']
strand_column = 6
start_column = 3
end_column = 4

extract_feature(records, species_table, output, id_column, feature_column, feature_type, strand_column,
            start_column, end_column, window_size)
