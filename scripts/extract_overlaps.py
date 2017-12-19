#!/usr/bin/env python3

"""This script will extract the pairwise overlapping nucleotides of each factors from the factor proxies
Computational time: ranges from 30 minutes to instant.
"""
from Bio import SeqIO
import sys
import glob
import os.path
from joblib import Parallel, delayed

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Species genome path:
species_genome = str(sys.argv[1])
# Species abbreviation:
species = '_'.join(str(species_genome.split('/')[-1]).split('_')[:2])
# Wanted number of threads at the same time:
n_threads = int(sys.argv[2])
# Follow up
follow_up = str(sys.argv[3])


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
# Count all the factor only record length
###
def count_pure_record_length(record, proxies_directory):
    record_ranges = proxies_directory + '/' + record.id

    length_factor_only = 0
    try:
        with open(record_ranges, 'r') as all_ranges:
            for each_range in all_ranges:
                line_range = each_range.strip().split()
                # Add the range length
                length_factor_only += len(range(int(line_range[0]), int(line_range[1]) + 1))
    # If no record for this factor, return length of 0
    except FileNotFoundError:
        return 0

    return length_factor_only


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
# Will check if we already computed this combination of factor/comparison but in a different order
###
def check_for_repeat(all_overlaps, factor, comparison):
    all_names = comparison.split('_')
    all_names.append(factor)
    repeat = False

    for each_already in all_overlaps:
        each_already = each_already.split('_')

        in_common = list(set(all_names) & set(each_already))

        if len(in_common) == len(all_names):
            repeat = True

    return repeat


def pairwise_overlap(record, factor_directory, factor, comparison_directory, comparison, overlap_directory):
    # Catch any missing record : as we will perform this step in parallel, we want to avoid raising errors
    if os.path.isfile(factor_directory + '/' + record.id) and os.path.isfile(comparison_directory + '/' + record.id):
        factor_ranges = extract_ranges(factor_directory + '/' + record.id)
        comparison_ranges = extract_ranges(comparison_directory + '/' + record.id)

        # Find all overlapping ranges
        overlap_ranges = intersects(factor_ranges, comparison_ranges)

        # If we do find some overlaps -> must write them to remove them from the factor ranges files later on
        if overlap_ranges:
            # We will now store the overlapping ranges
            checking_parent(overlap_directory + '/*')
            with open(overlap_directory + '/' + record.id, 'w') as overlap_file:
                for each_range in overlap_ranges:
                    overlap_file.write('\t'.join([str(each) for each in each_range]) + '\n')

            # We also return the factor_comparison to avoid repeats
            return '-'.join([factor, comparison])


def comparing_ranges(factor_to_compare, all_factors, records, species_proxies_directory):
    # We will store all the overlap for this pairwise comparisons
    all_overlaps = set()

    for each_factor in range(len(all_factors)):
        factor = all_factors[each_factor]
        for each_comparison in range(len(factor_to_compare)):
            comparison = factor_to_compare[each_comparison]

            # Avoid comparing factor overlap to himself...
            # Note the split by _ for further comparison (e.g. CDS_RNA vs CDS)
            if factor in comparison.split('_'):
                continue
            # We only do the "lower diagonal" when first comparing each individual factor to each other
            # as we don't want to do twice the same operation (e.g. CDS_RNA and RNA_CDS)
            elif check_for_repeat(all_overlaps, factor, comparison):
                continue
            else:
                factor_directory = species_proxies_directory + '/' + factor
                comparison_directory = species_proxies_directory + '/' + comparison
                overlap_directory = species_proxies_directory + '/' + '-'.join([factor, comparison])

                # Find the sum of all nucleotide which are either factor or comparison factor
                pure_factor_length = sum([count_pure_record_length(records[each_record], factor_directory)
                                          for each_record in range(len(records))])
                pure_comparison_length = sum([count_pure_record_length(records[each_record], comparison_directory)
                                              for each_record in range(len(records))])

                # Make sure neither are empty
                if pure_factor_length and pure_comparison_length:
                    overlaps = Parallel(n_jobs=n_threads)(delayed(pairwise_overlap)
                                                      (records[each_record], factor_directory, factor,
                                                       comparison_directory, comparison, overlap_directory)
                                                      for each_record in range(len(records)))
                    # If there did was an overlap in any of the records, keep it to avoid repeats
                    if any(overlaps):
                        all_overlaps.add([each for each in overlaps if each][0])

    return list(all_overlaps)


###
# Remove the overlapping ranges of factors
# Inputs:
#   - genomes_directory : directory of the genomes of species
#   - main_proxies_directory : main directory of the proxies obtained through the factor_proxies script
#   - all_factors : list of all wanted factors
#   - species : name of the wanted species
###
def extract_overlap(species_proxies_directory, all_factors, records):
    # This vector will help us track when we still have to make pairwise comparison of overlapping ranges
    to_compare = len(all_factors)

    # The first pairwise comparison we will make is all_factor VS all_factor anyway
    factor_to_compare = all_factors
    while to_compare >= 2:
        factor_to_compare = comparing_ranges(factor_to_compare, all_factors, records, species_proxies_directory)

        # If we have at least 2 set of overlapping factor, there is a chance nucleotide might be overlapping
        # in all of these different factors (e.g. if CDS_RNA and CDS_tandem, we might have CDS_RNA_tandem overlaps)
        to_compare = len(factor_to_compare)


# Fetch this species records for the lengths
records = fetch_fasta(species_genome)

# Directory which lead to the different species proxies
species_proxies_directory = '/'.join(["../files/factor_proxies", species])

# Find all the different factors ready
all_factors = [os.path.basename(each) for each in extract_path(species_proxies_directory + '/', '*')]

extract_overlap(species_proxies_directory, all_factors, records)

# Follow the progression of the analysis
checking_parent(follow_up)
with open(follow_up, 'w') as file:
    file.write('')
