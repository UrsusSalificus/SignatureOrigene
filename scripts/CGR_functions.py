# This script will compute the Chaos Game Representation (CGR) of a sequence.
from Bio import SeqIO
import sys
import math
import glob
import os
import numpy
from numpy import fft


###
# Check if parent directory is present, if not create it
###
def checking_parent (file_path):
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
# Extract all the path of files matching a certain pattern in a directory
###
def extract_path(files_directory, pattern):
    return glob.glob(files_directory + pattern)


###
# To makes up with any updates, we use the list of species from the downloading script:
###
def get_species(species_file='download_genomes.sh'):
    species=[]
    with open(species_file, 'r') as list_of_species:
        for each_line in list_of_species:
            if each_line.split('=')[0] == 'species':
                for each_species in each_line.split('=')[1].strip().strip('\'').split(' '):
                    species.append(each_species)
    return(species)


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
        checking_parent(outfile)
        with open(outfile, 'w') as file:
            for each_char in range(len(records)):
                file.write(str(xcord[each_char]) + '\t')
                file.write(str(ycord[each_char]) + '\n')

    # If no 2nd argument was given, outfile is empty (= considered False)
    else:
        # Store the CGR in a form of a list of list
        coordinates = [xcord, ycord]
        return coordinates


###
# Compute the k-mer Frequencies using the Chaos Game Representation (FCGR)
# Inputs:
#   - k_size : k-mer size
#   - CGR :
#       - Either a file (string)
#       - Or a set a of coordinates (list) obtained through CGR_coordinates function
#   - outfile : path to the output file, which will contain the FCGR
#           Note: if empty, will return the FCGR instead of writing a file.
# Output:
#   - Either a file, were each k-mer frequencies are separated by \t
#   - Or the k-mer frequencies stocked as a list
###
def FCGR_from_CGR(k_size, CGR, outfile):
    #####################
    # 1)  Fetch the coordinates and compute all the boundaries of the grid (for a certain k-mer size)
    #####################

    # If CGR is a string, it must be a path leading to a file containing all the coordinates
    if isinstance(CGR, str):
        x_coord = []
        y_coord = []
        with open(CGR, 'r') as CGR_file:
            for each_coord in CGR_file:
                x_coord.append(float(each_coord.split()[0]))
                y_coord.append(float(each_coord.split()[1]))
        coordinates = [x_coord, y_coord]
    # Else it's a Python list of coordinates
    else:
        coordinates = CGR
    # We take out the k_size-1 first coordinates, are these are the coordinates of words smaller than our
    # wanted k-mer size, and would add small errors later on when counting the frequencies
    for each_xy in [0, 1]:
        coordinates[each_xy] = coordinates[each_xy][k_size - 1::]

    # Now compute the different grid starting boundaries
    # Make sure the decimals are precise enough
    decimals = int(math.pow(10, k_size + 2))
    # Calculate the number of different k-mer we will compute (= number of grid we divide the CGR with)
    grid_size = int(math.pow(4, k_size))
    # Compute all the possible starting coordinates of all these grids
    start_coord_ranges = range(0, int(decimals), int(decimals / math.sqrt(grid_size)))

    #####################
    # 2)  After sorting the x coordinates, go one coordinates at a time, and check if it is a starting boundary
    #####################
    # For x coordinates:
    sort_x = sorted(coordinates[0])
    x_boundaries = []
    which_boundary = 0

    # We can go one sorted x coordinate at a time now:
    for each_coordinates in range(len(sort_x)):
        # If we are a the last coordinate, we must mark any remaining boundaries as empty, and we can stop here (break)
        if each_coordinates == len(sort_x) - 1:
            x_boundaries.append(each_coordinates)
            while len(x_boundaries) != len(start_coord_ranges) + 1:
                x_boundaries.append('empty')
            break
        # First, we see if the i coordinates is a "boundary" of a grid
        # (if bigger/outside of the i grid (represented by 'which_boundary'))
        if sort_x[each_coordinates] >= start_coord_ranges[which_boundary] / decimals:
            # If we are just before the last grid, we cannot have problem of empty grid left (see next if):
            if which_boundary == len(start_coord_ranges) - 1:
                x_boundaries.append(each_coordinates)
                which_boundary += 1
            # In some cases, the i coordinates is not only bigger than the i boundary, but also of i+1
            # (and even the next, so on and so on). In this case, we know that the i boundary
            # does not contain any coordinates ('empty'), and that we can go to the next boundary (which_boundary +1)
            elif sort_x[each_coordinates] >= start_coord_ranges[which_boundary + 1] / decimals:
                while not which_boundary == len(start_coord_ranges) - 1 \
                        and sort_x[each_coordinates] >= start_coord_ranges[which_boundary + 1] / decimals:
                    x_boundaries.append('empty')
                    which_boundary += 1
                x_boundaries.append(each_coordinates)
                which_boundary += 1
            # If not at the before-last grid, or if only bigger to i grid, it is the boundary of the i grid.
            # We then note where we are in the index, and go for the next starting boundary (which_boundary +1)
            else:
                x_boundaries.append(each_coordinates)
                which_boundary += 1
            # If we ended up being at the last grid, the last boundary is simply the last element
            # and we can stop here (break)
            if which_boundary == len(start_coord_ranges):
                x_boundaries.append(len(sort_x))
                break

    #####################
    # 3)  Quite similar step is done for each y coordinates corresponding to each x grids
    #####################
    # We will need the order of indices of x coordinates to do the same operations on y:
    numpy_x = numpy.array(coordinates[0])
    numpy_y = numpy.array(coordinates[1])
    sort_index_x = numpy.argsort(numpy_x)
    y_boundaries = []

    # Can go one starting boundary of grids of x at a time:
    # Note, we don't go for the -1 one, as it is not really a grid, but simply the ultimate grid's ending boundary
    for each_column in range(len(x_boundaries) - 1):
        y_boundaries.append([])

        # We will compute the y boundaries differently depending on whether the x column is empty or not:
        if not x_boundaries[each_column] == 'empty':
            # For each column of grids, we can sort the y of the sorted x, and do the same steps as before
            # First, we must find all the y corresponding to the x in the right grid:

            # In case we are at the last column
            if each_column == len(x_boundaries) - 1:
                corresponding_x_indexes = [x_boundaries[each_column], len(sort_x)]

            # Now in cases next column is an empty column, need the first non-empty one to know the corresponding x
            next_column = each_column + 1
            if x_boundaries[next_column] == 'empty':
                while True:
                    # In cases we have only empty columns left,
                    if next_column == len(x_boundaries) - 1:
                        corresponding_x_indexes = [x_boundaries[each_column], len(sort_x)]
                        break
                    # Else we must find which is the next non-empty
                    if x_boundaries[next_column] == 'empty':
                        next_column += 1
                    else:
                        corresponding_x_indexes = [x_boundaries[each_column], x_boundaries[next_column]]
                        break
            # Otherwise, it's simply the x between the boundaries
            else:
                corresponding_x_indexes = [x_boundaries[each_column], x_boundaries[next_column]]

            # Now that we have the corresponding x indexes, we can find the y of these corresponding x, sort them
            # and do the same operation that was performed on x before
            each_column_y = numpy_y[sort_index_x[corresponding_x_indexes[0]:corresponding_x_indexes[1]]]
            sort_y = numpy.sort(each_column_y)
            which_boundary = 0
            # Now do the same operation (see comments for the x coordinates operation)
            # Note, we do not place the last coordinate break at the beginning here,
            # as it would impede with counting later on.
            for each_coordinates in range(len(sort_y)):
                if sort_y[each_coordinates] >= start_coord_ranges[which_boundary] / decimals:
                    # Penultimate grid
                    if which_boundary == len(start_coord_ranges) - 1:
                        y_boundaries[each_column].append(each_coordinates)
                        which_boundary += 1
                    # Empty grid later on:
                    elif sort_y[each_coordinates] >= start_coord_ranges[which_boundary + 1] / decimals:
                        while not which_boundary == len(start_coord_ranges) - 1 \
                                and sort_y[each_coordinates] >= start_coord_ranges[which_boundary + 1] / decimals:
                            y_boundaries[each_column].append('empty')
                            which_boundary += 1
                        y_boundaries[each_column].append(each_coordinates)
                        which_boundary += 1
                    # Any not special grid
                    else:
                        y_boundaries[each_column].append(each_coordinates)
                        which_boundary += 1
                    # Ultimate grid:
                    if which_boundary == len(start_coord_ranges):
                        y_boundaries[each_column].append(len(sort_y))
                        break
                # If we are at the last coordinate
                if each_coordinates == len(sort_y) - 1:
                    y_boundaries[each_column].append(len(sort_y))
                    while len(y_boundaries[each_column]) != len(start_coord_ranges) + 1:
                        y_boundaries[each_column].append('empty')
                    break

        # Else, the whole column is empty, and must be marked as such:
        else:
            while len(y_boundaries[each_column]) != len(start_coord_ranges) + 1:
                y_boundaries[each_column].append('empty')

    #####################
    # 4)  Using only these y starting boundaries (indexes of y), we can count
    # how many coordinates there is in each boundary
    #####################
    FCGR = []
    # We know use the y boundaries to count the frequencies:
    for each_column in range(len(y_boundaries)):
        for each_kmer in range(len(y_boundaries[each_column]) - 1):
            # next_kmer will help us know which boundary we will use to compare to our actual kmer
            next_kmer = each_kmer + 1
            # If empty, this mean there is no counts for this kmer.
            if y_boundaries[each_column][each_kmer] == 'empty':
                FCGR.append(0)
            # As we subtract, we must have a non-empty index to compare with:
            elif y_boundaries[each_column][each_kmer + 1] == 'empty':
                find_non_empty = True
                while y_boundaries[each_column][next_kmer] == 'empty':
                    # If the ultimate kmer is also empty it simply means that each_kmer was the last boundary
                    if next_kmer == len(start_coord_ranges):
                        FCGR.append(0)
                        find_non_empty = False
                        break
                    # Else, we must keep on searching for the next which is not empty
                    else:
                        next_kmer += 1
                # If we find non-empty boundaries in the next kmers, can use it to count the number of kmer in grid i
                # Note that we add +1 to the index, as we are not as a starting position
                if find_non_empty:
                    # Note: Python start indexing at 0.
                    # To easily count using subtraction we must add +1 to the next index to get the real index.
                    FCGR.append((y_boundaries[each_column][next_kmer]) - y_boundaries[each_column][each_kmer])
            # Else can simply use the next kmer index to count the number of coordinates in the grid
            else:
                FCGR.append((y_boundaries[each_column][next_kmer]) - y_boundaries[each_column][each_kmer])

    # If outfile is non-empty, write the output
    if outfile:
        checking_parent(outfile)
        with open(outfile, 'w') as file:
            for each_count in FCGR:
                file.write(str(each_count) + '\t')
            file.write('\n')
        return FCGR
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


###
# Compute the power spectrum of the Chaos Game Representation (CGR) of a sequence,
# through the Discrete Fourier Transform (DFT) of the CGR
# Inputs:
#   - CGR :
#       - Either a file (string)
#       - Or a set a of coordinates (list) obtained through CGR_coordinates function
#   - outfile : path to the output file, which will contain the power spectrum
#           Note: if empty, will return the power spectrum instead of writing a file.
# Output:
#   - Either a file, were each power spectrum samples are separated by \t
#   - Or the power spectrum samples stocked as a list
###
def DFT_from_CGR(CGR, outfile):
    # If CGR is a string, it must be a path leading to a file containing all the coordinates
    if isinstance(CGR, str):
        x_coord = []
        y_coord = []
        with open(CGR, 'r') as CGR_file:
            for each_coord in CGR_file:
                x_coord.append(float(each_coord.split()[0]))
                y_coord.append(float(each_coord.split()[1]))
        coordinates = [x_coord, y_coord]
    # Else it's a Python list of coordinates
    else:
        coordinates = CGR

    # Z coordinates, as complex number of x and y coordinates
    z_coord = []
    for each_coord in range(len(coordinates[0])):
        z_coord.append(complex(coordinates[0][each_coord], coordinates[1][each_coord]))

    # Discrete Fourier Transform (DFT) through Fast Fourier Transform (FFT)
    # Will linearly decompose the CGR coordinates (the signal) into its component frequencies
    z_dft = fft.fft(z_coord)

    # Now we can compute the Power spectrum of the DFT = 'which frequencies contain the signal´s power'
    z_ps = []
    for each_coord in range(len(z_dft)):
        z_ps.append(abs(z_dft[each_coord]) ** 2)

    # If outfile is non-empty, write the output
    if outfile:
        checking_parent(outfile)
        with open(outfile, 'w') as file:
            for each_sample in z_ps:
                file.write(str(each_sample) + '\t')
            file.write('\n')
        return z_ps
    # If no 2nd argument was given, outfile is empty (= considered False)
    else:
        return z_ps