#!/usr/bin/env python3

"""This script compute the Pearson correlation distance between sets of FCGRs
"""
import glob
import sys
import os
import numpy
from scipy import stats as ss

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Path to the file containing the genomic signatures from all the regions concatenated:
concatenate = str(sys.argv[1])
# Species genome path:
species_genome = str(sys.argv[2])
# Species abbreviation:
species = '_'.join(str(species_genome.split('/')[-1]).split('_')[:2])
# Wanted window size:
window_size = int(sys.argv[3])
# Output file:
output = str(sys.argv[4])


###
# Extract all the path of files matching a certain pattern in a directory
###
def extract_path(files_directory, pattern):
    return glob.glob(files_directory + pattern)


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
# Compute the pairwise Pearson distances of a set of genomic signatures from different regions
# Inputs:
#   - genomic_signatures : path to a file containing all the region's genomic signatures
#   - n_region : number of region (lines in the genomic_signatures file)
#   - distance_matrix : path to the output file, which will contain the euclidean distance matrix
#           Note: if empty, will return the euclidean distance matrix instead of writing a file.
# Output:
#   - Either a file, where each distance between two regions are separated by \t,
#       forming a matrix of line separated by \n
#   - Or the Pearson distance matrix stocked as a numpy.array
# Performance :
#   - FCGRs : takes around 1 second for a set of 150'000 FCGRs (= 1 region).
###
def pearson_distance(genomic_signatures, n_region, distance_matrix):
    distances = []
    # We will have to compare line by line, in a pairwise fashion
    with open(genomic_signatures, 'r') as FCGRs:
        # i will help us know which line we have passed
        i = 0
        for each_region in FCGRs:
            # each_region = each line in the document
            distances.append([])
            # Note that when importing each line, we must exclude the first element (record id) and the last (empty)
            actual_line = numpy.array([float(each_DFT) for each_DFT in each_region.split('\t')[1:-1]])
            pair_line = numpy.array([float(each_DFT) for each_DFT in each_region.split('\t')[1:-1]])
            # j will help us know which line still have to do
            j = 0
            with open(genomic_signatures, 'r') as pair_FCGRs:
                # While we are at a pair already computed before, continue reading until find pair didn't compared
                while not j == i + 1:
                    # This will basically fill the lower diagonal of the distance matrix with 0
                    distances[i].append(0)
                    # Continue reading:
                    j += 1
                    pair_FCGRs.readline()
                # For each region left, compute the distance:
                for each_left in range(n_region - j):
                    # We import the line j to compare with i
                    pair_line = numpy.array([float(each_DFT) for each_DFT in pair_FCGRs.readline().split('\t')[1:-1]])
                    # Actual Pearson correlation distance computation (scipy.stats.pearsonr):
                    distances[i].append(ss.pearsonr(actual_line, pair_line)[0])     # Only take the r value [0]
            # This region/line done, update i
            i += 1

    # If outfile is non-empty, write the output
    if distance_matrix:
        checking_parent(distance_matrix)
        with open(distance_matrix, 'w') as outfile:
            for each_row in distances:
                for each_column in each_row:
                    outfile.write(str(each_column) + '\t')
                outfile.write('\n')
        return distances
    # If no 3rd argument was given, outfile is empty (= considered False)
    else:
        return distances


checking_parent(output)
# Total number of region (thus sum of all regions, in all records) :
n_region = sum([len(extract_path(each_record + '/', '*')) for each_record in extract_path(
                                '/'.join(['files/CGRs', str(window_size), species]) + '/', '*')])

pearson_distance(concatenate, n_region, output)

