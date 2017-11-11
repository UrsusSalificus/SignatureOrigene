#!/usr/bin/env python3

""" Extract introns from the gff file of the species
"""
import sys
import os

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# gff file path of the species:
species_table = str(sys.argv[1])
# Output path:
output = str(sys.argv[2])


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
# Add the introns in between exons in a gff file
# Inputs:
#   - species_table : path to the gff file of features (tested only on NCBI assemblies gff files)
#   - *_column : column number when separating feature line by \t
#   - outfile : path to the output file
# Output: gff file + the introns
###
def adding_introns(species_table, output, feature_column, strand_column,
                start_column, end_column):
    # Checking parent directory of output are present
    checking_parent(output)

    with open(species_table, 'r') as feature_table, open(output, 'w') as outfile:
        for each_line in feature_table:
            actual_line = each_line.split('\t')
            if len(actual_line) != 9:
                # We will use the fact that there is no \t in the comments to detect them
                outfile.write(each_line)
                # We store this line for the next line as a phony variable as it is not with the right structure
                previous_line = ['.'] * 9
            elif actual_line[feature_column] != 'exon' or actual_line[strand_column] != '+':
                # We also avoid all the lines which are not exons, or on the negative strand
                outfile.write(each_line)
                # Here, no phony, as it is already at the right structure
                previous_line = actual_line
            else:
                # We have stored the previous line each time, we can thus check whether it is an exon or not
                if previous_line[feature_column] == 'exon':
                    # We are in between two exons = intron! We will create the line to write:
                    start = int(previous_line[end_column]) + 1 # +1 to take intron nucleotide only
                    end = int(actual_line[start_column]) - 1
                    # First 2 elements and last 4 elements are the same for exon/intron
                    line_to_write = actual_line[0:2] + ['intron', start, end] + actual_line[5:]
                    for each_element in range(len(line_to_write)):
                        # Avoid the \t at the last element
                        outfile.write(str(line_to_write[each_element]) + '\t') if each_element != 8 \
                            else outfile.write(str(line_to_write[each_element]))
                        # Note: the last element always contain a \n at the end
                    # We also need the exon line written, and the previous line stored
                    outfile.write(each_line)
                    previous_line = actual_line
                else:
                    # If it is not, we are at the first exon -> just move on
                    outfile.write(each_line)
                    previous_line = actual_line


# Particularities of the feature table (e.g. which column contains what information):
feature_column = 2
# Indeed, introns are not directly annotated, but can be infered through exons
feature_type = 'exon'
strand_column = 6
start_column = 3
end_column = 4

adding_introns(species_table, output, feature_column, strand_column, start_column, end_column)
