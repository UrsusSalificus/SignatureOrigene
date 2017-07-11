# This script compute CGR of all windows of a certain size, in all species
from Bio import SeqIO
import glob
import sys
import os
from numpy import fft
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


for each_species in range(len(species)):
    DFT_directory = '/'.join(['../files/DFTs', str(window_in_kb), species[each_species]])
    concatenated_FCGRS = DFT_directory + '_DFTs'
    checking_parent(concatenated_FCGRS)
    with open(concatenated_FCGRS, 'w') as outfile:
        CGR_directory = '/'.join(['../files/CGRs', window_in_kb, species[each_species]]) # Path to directory
        all_records = extract_path(CGR_directory + '/', '*')
        for each_record in range(len(all_records)):
            CGR_files = extract_path(str(all_records[each_record] + '/'), '*')
            record_name = all_records[each_record].split('/')[-1]
            DFT_region = '/'.join([DFT_directory, record_name, 'DFT_region_'])   # Path to region
            DFTs = Parallel(n_jobs=n_threads)(delayed(DFT_from_CGR)
                                               (CGR_files[each_region], DFT_region + str(each_region))
                                               for each_region in range(len(CGR_files)))
            for each_region in range(len(DFTs)):
                outfile.write(record_name + '\t')
                for each_count in DFTs[each_region]:
                    outfile.write(str(each_count) + '\t')
                outfile.write('\n')





