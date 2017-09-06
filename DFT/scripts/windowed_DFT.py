#!/usr/bin/env python3

"""This script compute k-mer frequencies using CGR of all windows of a certain size
"""
import glob
import sys
import os
import numpy
from joblib import Parallel, delayed

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Species abbreviation:
species = str(sys.argv[2])
# Wanted window size:
window_size = int(str(sys.argv[3]).split('/')[-1])
# Wanted number of threads at the same time:
n_threads = int(sys.argv[4])
# Output file:
output = str(sys.argv[5])


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


checking_parent(output)
# Opening concatenated file on top level, to avoid rewriting at each record
with open(output, 'w') as outfile:
    # Path to CGR directory
    CGR_directory = '/'.join(['../files/CGRs/', str(window_size), species])
    # Get all the different CGRs files path
    all_records = extract_path(CGR_directory + '/', '*')
    for each_record in range(len(all_records)):
        # Find all CGR files with glob:
        CGR_files = extract_path(str(all_records[each_record] + '/'), '*')
        # Extract all the record names and store it for later:
        record_name = all_records[each_record].split('/')[-1]
        # Parallel computation for every region, stored in individual files:
        DFTs = Parallel(n_jobs=n_threads)(delayed(DFT_from_CGR)(CGR_files[each_region], '')
                                          for each_region in range(len(CGR_files)))

        # Write each region's genomic signature in a single concatenated file:
        for each_region in range(len(DFTs)):
            outfile.write(record_name + '\t')
            for each_count in DFTs[each_region]:
                outfile.write(str(each_count) + '\t')
            outfile.write('\n')
