#!/usr/bin/env python3

"""This script contain function built but not used directly in the analysis
"""
from Bio import SeqIO
import sys
import math
import glob
import os
import numpy
from numpy import fft

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"


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
# Compute the indexes of the k-mer Frequencies of the Chaos Game Representation (FCGR)
# Inputs:
#   - k_size : k-mer size
#   - halving: whether or not (by default = False) we want the indexes of the sum of the reverse complementary kmer
# Output:
#   - The k-mer indexes (length(indexes) = all possible k-mer)
###
def FCGR_indexes(k_size, halving = False):
    # Get the indexes for this specific k-mer size
    # Will take advantage of the fact that dividing the picture in grid yield mini tetramer square of indexes:
    # Basically, a grid in an odd column and odd line yields a A (0,0), odd column and even line yields a C (0,1),
    # even column and odd line yields a T (1,0) and even column and even line yields a G (1,1)
    # First, we see what is the maximum number of mini-square we have in the picture
    grid_size = int(math.pow(4, k_size))
    mini_square = int(math.sqrt(grid_size))
    grid_indexes = []
    # For each k-mer size, find the indexes
    for each_size in range(k_size):
        # As python range begin at 0, must +1 to get to real k-mer size
        each_size += 1
        # actual_grid_index will contain all the indexes for this k-mer size
        actual_grid_index = [''] * grid_size
        # This pointer will track movement on the grid
        grid_pointer = 0
        # This will give us the number of grid for this size,
        # which is used to find the number of mini-squares for this size
        actual_grid_size = int(math.pow(4, each_size))
        actual_mini_square = int(math.sqrt(actual_grid_size))
        # multiplier will help use reach the max ammount of grid whenever we are at < size than max k-mer size
        # Note: multiplying by 1 whenever we are at same size than max k-mer size
        multiplier = mini_square / actual_mini_square
        # We have to tweak the column/lines to make it works whenever we are at < size than max k-mer size
        actual_lines = []
        actual_column = []
        # This part will thus only change line/column whenever we are at < size than max k-mer size
        for each_mini_square in range(actual_mini_square):
            multiplied_line = [each_mini_square] * int(multiplier)
            multiplied_column = [each_mini_square] * int(multiplier)
            for each_new_line in multiplied_line:
                actual_lines.append(each_new_line)
            for each_new_column in multiplied_column:
                actual_column.append(each_new_column)

        # Now can see whether odd/even column and stock the indexes for this size
        for each_column in actual_column:
            if not (each_column + 1) % 2 == 0:
                for each_line in actual_lines:
                    if not (each_line + 1) % 2 == 0:
                        actual_grid_index[grid_pointer] = 'A'
                        grid_pointer += 1
                    else:
                        actual_grid_index[grid_pointer] = 'C'
                        grid_pointer += 1
            else:
                for each_line in actual_lines:
                    if not (each_line + 1) % 2 == 0:
                        actual_grid_index[grid_pointer] = 'T'
                        grid_pointer += 1
                    else:
                        actual_grid_index[grid_pointer] = 'G'
                        grid_pointer += 1
        grid_indexes.append(actual_grid_index)

    # Lastly, we now have k lists of indexes, which we must join to have the real index per grid
    indexes = []
    # If k = 1, we don't have to go list by list (only one list of index)
    if len(grid_indexes) > 1:
        for each_index_grid in range(len(grid_indexes[0])):
            single_index = []
            for each_index_size in reversed(range(len(grid_indexes))):
                single_index.append(grid_indexes[each_index_size][each_index_grid])
            indexes.append(''.join(single_index))
    else:
        indexes = grid_indexes[0]

    if halving:
        # To have reverse complementary kmer joined, we must kind of cut the vector in 2, then join the right column
        # to one another (first to last, second to last-1, and so on)
        cut_indexes = list()
        for each_half_column in range(int(mini_square / 2)):
            positive_strand_start = each_half_column * mini_square
            positive_strand_end = positive_strand_start + mini_square
            positive_indexes = indexes[positive_strand_start:positive_strand_end]
            negative_strand_start = ((mini_square - each_half_column) - 1) * mini_square
            negative_strand_end = negative_strand_start + mini_square
            negative_indexes = indexes[negative_strand_start:negative_strand_end]

            cut_indexes.extend(['/'.join(x) for x in zip(*[positive_indexes, negative_indexes])])

        indexes = cut_indexes

    return indexes


###
# Compute the pairwise euclidean distances of a set of genomic signatures from different regions
# Inputs:
#   - genomic_signatures : path to a file containing all the region's genomic signatures
#   - n_region : number of region (lines in the genomic_signatures file)
#   - distance_matrix : path to the output file, which will contain the euclidean distance matrix
#           Note: if empty, will return the euclidean distance matrix instead of writing a file.
# Output:
#   - Either a file, where each distance between two regions are separated by \t,
#       forming a matrix of line separated by \n
#   - Or the euclidean distance matrix stocked as a numpy.array
# Performance :
#   - DFTs : takes around 1 second for a set of 150'000 DFTs (= 1 region).
###
def euclidean_distance(genomic_signatures, n_region, distance_matrix):
    distances = []
    with open(genomic_signatures, 'r') as DFTs:
        i = 0
        for each_region in DFTs:
            distances.append([])
            # Note that when importing each line, we must exclude the first element (record id) and the last (empty)
            actual_line = numpy.array([float(each_DFT) for each_DFT in each_region.split('\t')[1:-1]])
            j = 0
            with open(genomic_signatures, 'r') as pair_DFTs:
                # While we are at a pair already computed before, continue reading until find pair didn't compared
                while not j == i + 1:
                    distances[i].append(0)
                    j += 1
                    pair_DFTs.readline()
                for each_left in range(n_region - j):
                    pair_line = numpy.array([float(each_DFT) for each_DFT in pair_DFTs.readline().split('\t')[1:-1]])
                    distances[i].append(numpy.linalg.norm(actual_line - pair_line))
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