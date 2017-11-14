#!/usr/bin/env python3

""" Masking to only have nucleotides composed of the wanted feature
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os
import numpy as np

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
# Output:
#   - A numpy array of n (number of different records in records) proxies
#       - Each proxy contain m (number of nucleotides in the i record) values
#       - Each value is either 0 to indicate that this factor is not involved in the wanted factor
#           or 1 to indicate that this factor is linked to the wanted factor
###
def extract_factor(records, factor_type, feature_type, species_table, id_column, feature_column,
                   strand_column, start_column, end_column):
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

            # Translate it to numpy array
            record_proxy = np.asarray(record_proxy)

            # Add the proxy of factor to the list
            proxies_records.append(record_proxy)
    return proxies_records


###
# Extract from the proxy all the nucleotide that are from the wanted factor
# Inputs:
#   - records : fetched sequence (fasta) of the species whole genome
#   - proxies_records : list of proxies obtained through the extract_factor function
#   - factor : name/abbreviation of the wanted factor
# Output:
#   - A list of SeqRecords of pure factor
###
def extract_pure_ranges(records, proxies_records, factor):
    # We will store the pure nucleotides records
    factor_records = list()

    for each_proxy in range(len(proxies_records)):
        # We will find all non-overlapping ranges of the pure factors:
        # Only difference with masked -> == 1
        no_overlap_ranges = np.where(proxies_records[each_proxy] == 1)

        # Translate the record in numpy array
        record_array = np.array(list(str(records[each_proxy].seq)))
        factor_only = ''.join(record_array[no_overlap_ranges])

        new_record = SeqRecord(seq=factor_only, id=factor)
        factor_records.append(new_record)

    return factor_records


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

# Compute factor proxy of records
proxies = extract_factor(records, factor_type, feature_type, species_table, id_column, feature_column,
                         strand_column, start_column, end_column)

# Compute records composed of only the factor's nucleotides
pure_records = extract_pure_ranges(records, proxies, factor)

# Checking parent directory of output are present
checking_parent(output)

# Writing the fasta file
with open(output, 'w') as outfile:
    # We might end up with records which are too small
    # In this case, we must concatenate them:
    if any([len(each_record) < window_size for each_record in pure_records]):
        concatenated_seq = str()
        for each_record in pure_records:
            concatenated_seq += each_record
        concatenated_records = SeqRecord(seq=concatenated_seq.seq, id=factor)

        # Write the new list of records
        SeqIO.write(concatenated_records, outfile, "fasta")
    else:
        # No need to concatenate -> directly write the new list of records
        SeqIO.write(pure_records, outfile, "fasta")
