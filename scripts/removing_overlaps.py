#!/usr/bin/env python3

"""This script will remove the overlapping nucleotides of two factors from the factor ranges
Computational time: ranges from 30 minutes (H. sapiens/M. musculus) to instant.
"""
from Bio import SeqIO
import sys
import glob
import os.path
import shutil
from joblib import Parallel, delayed

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Species genome path:
species_genome = str(sys.argv[1])
# Species abbreviation:
species = '_'.join(str(species_genome.split('/')[-1]).split('_')[:2])
# Wanted number of threads at the same time:
n_threads = int(sys.argv[3])
# Follow up
follow_up = str(sys.argv[4])


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
# Extract ranges as a list of ranges
###
def extract_ranges(ranges_file_path):
    ranges = list()
    with open(ranges_file_path, 'r') as ranges_file:
        for each_line in ranges_file:
            ranges.append([int(each_line.strip().split('\t')[0]), int(each_line.strip().split('\t')[1])])
    return ranges


###
# Find intersects between two list of ranges (a and b)
###
def intersects(a, b):
    overlap_ranges = list()
    i = j = 0
    while i < len(a) and j < len(b):
        a_left, a_right = a[i]
        b_left, b_right = b[j]
        end_pts = sorted([a_left, a_right, b_left, b_right])
        middle = [end_pts[1], end_pts[2]]

        if a_right < b_right:
            i += 1
        else:
            j += 1

        if a_right >= b_left and b_right >= a_left:
            overlap_ranges.append(middle)

    length_of_list = len(overlap_ranges)

    # We will also merge two intervals which are next to each other
    i = 0
    while i < (length_of_list - 1):
        if overlap_ranges[i][1] in (overlap_ranges[i + 1][0], overlap_ranges[i + 1][0] - 1):
            overlap_ranges[i:i + 2] = [[overlap_ranges[i][0], overlap_ranges[i + 1][1]]]
            length_of_list -= 1
        else:
            i += 1

    return overlap_ranges


###
# Clean factor ranges of overlapping ranges
###
def cleaning_ranges(record, factor_directory, overlap_directory):
    # Catch any missing record : as we will perform this step in parallel, we want to avoid raising errors
    if os.path.isfile(factor_directory + '/' + record.id) and os.path.isfile(overlap_directory + '/' + record.id):
        factor_ranges = extract_ranges(factor_directory + '/' + record.id)
        overlap_ranges = extract_ranges(overlap_directory + '/' + record.id)

        # As the removing is a multilevel process, we may have some overlap which may have been already erased
        # In this case, we must find which are the overlap to still remove
        overlap_ranges = intersects(factor_ranges, overlap_ranges)

        i = j = 0
        # Up to the last overlapping range
        while j < len(overlap_ranges):
            a_left, a_right = factor_ranges[i]
            b_left, b_right = overlap_ranges[j]

            end_pts = sorted([a_left, a_right, b_left, b_right])
            bottom = [end_pts[0], end_pts[1]]
            tail = [end_pts[2], end_pts[3]]

            if a_right < b_right:
                i += 1
            else:
                j += 1

            if a_right >= b_left and b_right >= a_left:
                # If everything is erased = whole range is gone:
                if bottom[0] == bottom[1] and tail[0] == tail[1]:
                    del factor_ranges[i]
                # If left some part at the bottom:
                elif bottom[0] != bottom[1] and tail[0] == tail[1]:
                    factor_ranges[i] = [bottom[0], bottom[1] - 1]
                # Or at the end:
                elif bottom[0] == bottom[1] and tail[0] != tail[1]:
                    factor_ranges[i] = [tail[0] + 1, tail[1]]
                # Or finally both:
                else:
                    factor_ranges[i] = [bottom[0], bottom[1] - 1]
                    factor_ranges.insert(i + 1, [tail[0] + 1, tail[1]])

        # If we still find some overlaps -> can write them
        if factor_ranges:
            # We will now store the clean ranges
            with open(factor_directory + '/' + record.id, 'w') as clean_file:
                for each_range in factor_ranges:
                    clean_file.write('\t'.join([str(each) for each in each_range]) + '\n')
        # If we lost all ranges, we will simply delete the directory
        else:
            if os.path.isdir(factor_directory):
                shutil.rmtree(factor_directory)


###
# Remove the overlapping ranges of factors
# Inputs:
#   - genomes_directory : directory of the genomes of species
#   - main_proxies_directory : main directory of the proxies obtained through the factor_proxies script
#   - all_factors : list of all wanted factors
#   - species : name of the wanted species
###
def remove_overlap(to_clean, species_proxies_directory):
    # We will perform this for each level (up to *number_of_individual_factors* - 1)
    for each_level in range(len(to_clean) - 1):
        factor_to_clean = to_clean[each_level]
        level_to_clean = to_clean[each_level + 1]

        # If there is only one element in clean, must make it a list
        if isinstance(level_to_clean, str):
            level_to_clean = [level_to_clean]

        for each_factor in range(len(factor_to_clean)):
            factor = factor_to_clean[each_factor]

            for each_overlap in range(len(level_to_clean)):
                overlap = level_to_clean[each_overlap]
                # We only remove when factor in the overlapping factors
                if len(list(set(factor.split('-')) & set(overlap.split('-')))) == len(factor.split('-')):
                    factor_directory = species_proxies_directory + '/' + factor
                    overlap_directory = species_proxies_directory + '/' + overlap

                    Parallel(n_jobs=n_threads)(delayed(cleaning_ranges)
                                               (records[each_record], factor_directory, overlap_directory)
                                               for each_record in range(len(records)))


# Fetch this species records for the lengths
records = fetch_fasta(species_genome)

# Directory which lead to the different species proxies
species_proxies_directory = '/'.join(["../files/factor_proxies", species])

# Find all the different factors ready
all_levels = [os.path.basename(each) for each in extract_path(species_proxies_directory + '/', '*')]

# How many different levels we have?
n_levels = max([len(each.split('-')) for each in all_levels])

to_clean = list()
# Note: we add + 1 due to the pythonic counting starting at 0
for each_level in range(1, n_levels + 1):
    to_clean.append([each for each in all_levels if len(each.split('-')) == each_level])

remove_overlap(to_clean, species_proxies_directory)

# Follow the progression of the analysis
checking_parent(follow_up)
with open(follow_up, 'w') as file:
    file.write('')
