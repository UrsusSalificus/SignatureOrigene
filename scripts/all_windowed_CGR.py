# This script compute CGR of all windows of a certain size, in all species
from Bio import SeqIO
import glob
import math
import sys
import os
from joblib import Parallel, delayed

# Wanted window size:
window_size = int(sys.argv[1])
window_in_kb = str(window_size)[:-3] + 'kb'
# Wanted number of threads at the same time:
n_threads = int(sys.argv[2])


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
# Slight update to make the computing CGR function to ignore sequences containing unknown nucleotides
###
def N_sensitive_CGR(window, outfile):
    # If any character in the sequence is NOT a standard nucleotides (including unknown nucleotides), do NOT compute:
    if not any([c not in 'ATCGatcg' for c in window]):
        CGR_coordinates(window, outfile)


# For each species, we will compute CGR on each window, of each records.
for each_species in range(len(species)):
    # Due to non-consistent pattern of file name, whole genome ('genomic') is used in multiple
    # fasta files (CDS or RNA only, and any nucleotides) names. We must thus reconstruct the exact path.
    # To do so, we will use the feature_table (only one per species)
    species_table = extract_path('../data/genomes/' + species[each_species], '*_feature_table*')[0]

    split_index = species_table.split('/')[3].split('_')
    # The index word length varies between species : we must skip the two first  word (species name)
    index = ''
    walking = 2
    while split_index[walking] != 'feature':
        index += split_index[walking] + '_'
        walking +=1

    pattern_genome = str(species[each_species] + '_' + index + 'genomic*')
    species_genome = extract_path('../data/genomes/', pattern_genome)[0]

    # Fetching the genomic fasta file
    records = fetch_fasta(species_genome)
    # Will work for single record too:
    for each_sequence in records:
        if len(each_sequence.seq) > window_size:
            n_windows = math.floor(len(each_sequence.seq) / window_size)  # Number of windows
            seq_directory = '/'.join(['../files/CGRs', window_in_kb, species[each_species],
                                      each_sequence.id, 'CGR_region_'])
            # Parallel the CGR on n_jobs core
            Parallel(n_jobs=n_threads)(delayed(N_sensitive_CGR)
                                       (str(each_sequence.seq[
                                            (start * window_size):((start * window_size) + window_size)]),
                                        seq_directory + str(start))
                                       for start in range(0, n_windows))

