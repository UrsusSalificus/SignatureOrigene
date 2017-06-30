# This script will compute the Chaos Game Representation (CGR) of a sequence.
from Bio import SeqIO
import sys
import math
import glob
import os


###
# Check if parent directory is present, if not create it
###
def checking_parent (file_path):
    # We don't want need the file name, so will take everything but the last part
    parent_directories = '/'.join(file_path.split('/')[0:(len(file_path.split('/')) - 1)])
    if not os.path.exists(parent_directories):
        os.makedirs(parent_directories)


###
# Extract all the path of files matching a certain pattern in a directory
###
def extract_path(files_directory, pattern):
    return glob.glob(files_directory + pattern)

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
    return(records)


###
# Compute the Chaos Game Representation (CGR) of the cleaned sequence
# Inputs:
#   - records : fetched sequence (fasta) one wants the CGR computed on
#   - outfile : path to the output file, which will contain the x/y CGR coordinates
#           Note: if empty, will return the coordinates instead of writing a file.
# Output:
#   - Either a file, were each line contain a set of x/y coordinates (separated by \t)
#   - Or the coordinates stocked as [[x coordinates], [y coordinates]]
###
def CGR_coordinates(records, outfile):
    # Prepare the list of x and y coordinates, already at the right size as we will us numeric range in our loop
    xcord = [0] * len(records)
    ycord = [0] * len(records)
    # actual variables are the pointer variables (were are we?), it will change at each nucleotide.
    actual_x = 0.5
    actual_y = 0.5
    # For each nucleotide, from your actual position in the picture, go half-way to this nucleotide corner
    # Nucleotide corner (x axis, y axis): A (0,0), T (1,0), C (0,1) and G (1,1)
    for each_char in range(len(records)):
        if records[each_char] == "A" or records[each_char] == "a":
            xcord[each_char] = actual_x + ((-actual_x) / 2)
            ycord[each_char] = actual_y + ((-actual_y) / 2)
        elif records[each_char] == "C" or records[each_char] == "c":
            xcord[each_char] = actual_x + ((-actual_x) / 2)
            ycord[each_char] = actual_y + ((1 - actual_y) / 2)
        elif records[each_char] == "G" or records[each_char] == "g":
            xcord[each_char] = actual_x + ((1 - actual_x) / 2)
            ycord[each_char] = actual_y + ((1 - actual_y) / 2)
        elif records[each_char] == "T" or records[each_char] == "t":
            xcord[each_char] = actual_x + ((1 - actual_x) / 2)
            ycord[each_char] = actual_y + ((-actual_y) / 2)
        else:
            raise ValueError('Nucleotide not valid (not ATCG), may impede the rest of the analysis, please clean the '
                             'sequence')
        actual_x = xcord[each_char]
        actual_y = ycord[each_char]

    # If outfile is non-empty, write the output
    if outfile:
        with open(outfile, 'w') as file:
            for each_char in range(len(records)):
                file.write(str(xcord[each_char])+ '\t')
                file.write(str(ycord[each_char]) + '\n')

    # If no 2nd argument was given, outfile is empty (= considered False)
    else:
        # Store the CGR in a form of a list of list
        coordinates = [xcord, ycord]

    return coordinates



###
# Compute the k-mer Frequencies of the Chaos Game Representation (FCGR) of the cleaned sequence
# Inputs:
#   - k_size : k-mer size
#   - CGR :
#       - Either a file (string)
#       - Or a set a of coordinates (list) obtained through CGR_coordinates function
#   - outfile : path to the output file, which will contain the FCGR
#           Note: if empty, will return the FCGR instead of writing a file.
# Output:
#   - Either a file, were each k-mer frequencies are separated by \t
#   - Or the coordinates stocked in the same format
###
def FCGR_from_CGR(k_size, CGR, outfile):
    # If CGR is a string = a path
    if isinstance(CGR, str):
        x_coord = []
        y_coord = []
        with open(CGR, 'r') as CGR_file:
            for each_coord in CGR_file:
                x_coord.append(each_coord.split()[0])
                y_coord.append(each_coord.split()[1])
        CGR = [x_coord, y_coord]
    # Else it can be used as it is

    # Make sure the decimals are precise enough
    decimals = int(math.pow(10, k_size + 2))
    # Calculate the number of different k-mer we will compute (= number of grid we divide the CGR with)
    grid_size = int(math.pow(4, k_size))

    # Compute all the possible starting/ending coordinates of all these grids
    start_coord_ranges = range(0, int(decimals), int(decimals / math.sqrt(grid_size)))
    end_coord_ranges = range(int(decimals / math.sqrt(grid_size)), int(decimals + decimals / math.sqrt(grid_size)),
                             int(decimals / math.sqrt(grid_size)))

    # Compute all the grid x/y coordinates under the format [x minimum, x maximum, y minimum, y maximum]
    # Note: grid_ranges have same length than grid_size
    grid_ranges = []
    for each_column in range(int(math.sqrt(grid_size))):
        # Each "column" of grid has same max/min x coordinates, but all the possible different max/min y coordinates
        x_start = start_coord_ranges[each_column] / decimals
        x_end = end_coord_ranges[each_column] / decimals
        for each_line in range(int(math.sqrt(grid_size))):
            y_start = start_coord_ranges[each_line] / decimals
            y_end = end_coord_ranges[each_line] / decimals
            grid_ranges.append([x_start, x_end, y_start, y_end])

    FCGR = []
    coordinates = list(CGR)
    # Compute the frequency of k-mer:
    #   by dividing the Chaos Game Representation (CGR) of the sequence into multiple grids
    # Note: we will intentionally exclude the k-1 first elements (as they are smaller than the k-mer size we want)
    for each_grid in grid_ranges:
        # grid_count will count the number of "points" in the actual grid
        grid_count = 0
        # We also store the coordinates that were already assigned, to reduce computational time
        to_remove = []
        # For each set of coordinates, check if in the grid
        for each_set in range(len(coordinates[0])):
            # If this set of coordinates is in the actual grid,
            # increase actual grid count by 1 and stock it to remove it
            if each_grid[0] < float(coordinates[0][each_set]) < each_grid[1] \
                    and each_grid[2] < float(coordinates[1][each_set]) < each_grid[3]:
                grid_count += 1
                to_remove.append(each_set)
        # Stock the output
        FCGR.append(grid_count)

        # Remove the assigned coordinates (to reduce computational time)
        # We sort the sets to remove to begin deleting by the end = no problem of moving indexes
        for each_rm in sorted(to_remove, reverse=True):
            del coordinates[0][each_rm]
            del coordinates[1][each_rm]

    # If outfile is non-empty, write the output
    if outfile:
        with open(outfile, 'w') as file:
            for each_char in coordinates:
                file.write(each_char + '\t')
            file.write('\n')

    # If no 2nd argument was given, outfile is empty (= considered False)
    else:
        return FCGR


###
# Compute the indexes of the k-mer Frequencies of the Chaos Game Representation (FCGR)
# Inputs:
#   - k_size : k-mer size
# Output:
#   - The k-mer indexes (length(indexes) = all possible k-mer)
###
def FCGR_indexes(k_size):
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
        # This will give us the number of grid for this size, which is used to find the number of mini-squares for this size
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
    # If k = 1, we don't have ot go list by list (only one list of index)
    if len(grid_indexes) > 1:
        for each_index_grid in range(len(grid_indexes[0])):
            single_index = []
            for each_index_size in reversed(range(len(grid_indexes))):
                single_index.append(grid_indexes[each_index_size][each_index_grid])
            indexes.append(''.join(single_index))
    else:
        indexes = grid_indexes[0]

    return indexes
