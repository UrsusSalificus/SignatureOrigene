# This script compute the power spectrum of the Chaos Game Representation (CGR) of a sequence,
# through the Discrete Fourier Transform (DFT) of a CGR
import glob
import sys
import os
import numpy
from joblib import Parallel, delayed

# Wanted window size:
window_size = int(sys.argv[1])
window_in_kb = str(window_size)[:-3] + 'kb'
# Wanted number of threads at the same time:
n_threads = int(sys.argv[2])
# Which genomic signature to use?
gs_type = int(sys.argv[3])


###
# To makes up with any updates, we use the list of species from the downloading script:
###
def get_species(species_file='download_genomes.sh'):
    species = []
    with open(species_file, 'r') as list_of_species:
        for each_line in list_of_species:
            if each_line.split('=')[0] == 'species':
                for each_species in each_line.split('=')[1].strip().strip('\'').split(' '):
                    species.append(each_species)
    return (species)


# Get all the species abbreviation for the study
species = get_species()


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

Parallel(n_jobs=n_threads)(delayed(euclidean_distance)
                           ('/'.join(['../files', gs_type, window_in_kb,
                                      '_'.join([species[each_species], gs_type])]),
                            sum([len(extract_path(each_record + '/', '*')) for each_record in extract_path(
                                '/'.join(['../files/CGRs', window_in_kb, species[each_species]]) + '/', '*')]),
                            '/'.join(['../files/distances/euclidian', '_'.join([window_in_kb, gs_type]),
                                      '_'.join([species[each_species], 'euc_dist'])])
                            )
                           for each_species in range(len(species)))

