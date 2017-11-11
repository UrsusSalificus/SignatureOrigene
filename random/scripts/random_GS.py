"""This script will compute the CGR, then both genomic signatures (DFT/FCGR) of n random windows of each species
"""
from Bio import SeqIO
import random
import sys
import math
import numpy
import os
import glob
import os.path
from joblib import Parallel, delayed


__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Wanted number of random windows per species
sample_size = int(sys.argv[1])
# Wanted window size:
window_size = int(sys.argv[2])
# Wanted k-mer size:
k_size = int(sys.argv[3])
# Wanted number of threads at the same time:
n_threads = int(sys.argv[4])
# Output FCGR file:
output_FCGRs = "temp/FCGRs"
# Output DFT file:
output_DFTs = "temp/DFTs"


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
# Compute the Chaos Game Representation (CGR) of a sequence
# Inputs:
#   - records : fetched sequence (fasta) one wants the CGR computed on
#   - outfile : path to the output file, which will contain the x/y CGR coordinates
#           Note: if empty, will return the coordinates instead of writing a file.
# Output:
#   - Either a file, where each line contain a set of x/y coordinates (separated by \t)
#   - Or the coordinates stocked as [[x coordinates], [y coordinates]]
# Performance :
#   - Takes around 1 second for a 300'000 base pairs long sequence.
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
#   - Either a file, where each k-mer frequencies are separated by \t
#   - Or the k-mer frequencies stocked ias a list
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
# Compute the power spectrum of the Chaos Game Representation (CGR) of a sequence,
# through the Discrete Fourier Transform (DFT) of the CGR
# Inputs:
#   - CGR :
#       - Either a file (string)
#       - Or a set a of coordinates (list) obtained through CGR_coordinates function
#   - outfile : path to the output file, which will contain the power spectrum
#           Note: if empty, will return the power spectrum instead of writing a file.
# Output:
#   - Either a file, where each power spectrum samples are separated by \t
#   - Or the power spectrum samples stocked as a list
# Performance :
#   - Takes around 1 second for a set of 150'000 x/y coordinates.
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
    z_dft = numpy.fft.fft(z_coord)

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
    # If no 2nd argument was given, outfile is empty (= considered False)
    else:
        return z_ps

# Extract all the species compared
species = [os.path.basename(each) for each in glob.glob('config/species/*')]

all_FCGRs = []
all_DFTS = []
checking_parent(output_FCGRs)
checking_parent(output_DFTs)
with open(output_FCGRs, 'w') as outfile_FCGRs, open(output_DFTs, 'w') as outfile_DFTs:
    for each_species in species:
        species_genome = "../data/genomes/" + each_species + "_genomes.fna"

        # Fetching the genomic fasta file
        records = fetch_fasta(species_genome)

        all_possibilities = []
        # We want to have the same chance for any window of the genome to be present
        # Thus, we sample randomly in all the possible window, in all the different records (all_possibilities)
        for each_record in range(len(records)):
            if len(records[each_record].seq) >= window_size:
                number_window = math.floor(len(records[each_record].seq)/window_size)
                all_window = ['-'.join([str(each_record), str(each_window)]) for each_window in range(number_window)]
                all_possibilities.extend(all_window)

        # Check there is no unknown nucleotides in the window
        # Why now and not before?
        # Because by picking 20 window, we greatly reduce risk of finding unknown nucleotide in seq
        ready_to_go = False
        while not ready_to_go:
            candidates = random.sample(all_possibilities, sample_size)
            all_window=[]
            for each_candidate in candidates:
                record_number = int(each_candidate.split('-')[0])
                window_number = int(each_candidate.split('-')[1])
                start = window_size*window_number
                end = start  + window_size
                all_window.extend(records[record_number].seq[start:end])
            if not any([c not in 'ATCGatcg' for c in all_window]):
                ready_to_go = True

        CGRs = Parallel(n_jobs=n_threads)\
            (delayed(CGR_coordinates)
             (records[int(each_candidate.split('-')[0])].seq
              [(window_size * int(each_candidate.split('-')[1])):
              ((window_size * int(each_candidate.split('-')[1])) + window_size)],
              '') for each_candidate in candidates)
        FCGRs = Parallel(n_jobs=n_threads)(delayed(FCGR_from_CGR)(k_size, each_CGR,'') for each_CGR in CGRs)
        DFTs = Parallel(n_jobs=n_threads)(delayed(DFT_from_CGR)(each_CGR,'') for each_CGR in CGRs)

        # Write each region's genomic signature in a single concatenated file:
        for each_candidate in range(len(candidates)):
            outfile_FCGRs.write(each_species + '\t')
            outfile_DFTs.write(each_species + '\t')
            for each_count in FCGRs[each_candidate]:
                outfile_FCGRs.write(str(each_count) + '\t')
            for each_count in DFTs[each_candidate]:
                outfile_DFTs.write(str(each_count) + '\t')
            outfile_FCGRs.write('\n')
            outfile_DFTs.write('\n')

