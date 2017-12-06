"""This script will extract the overall percentages of each factor, in each species
"""
from Bio import SeqIO
import sys
import glob
import os.path
import math
from joblib import Parallel, delayed

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Wanted window size:
window_size = int(sys.argv[1])
# Wanted number of threads at the same time:
n_threads = int(sys.argv[2])
# Wanted window size:
output = str(sys.argv[3])


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
    with open(record_ranges, 'r') as all_ranges:
        for each_range in all_ranges:
            line_range = each_range.strip().split()
            # Add the range length
            length_factor_only += len(range(int(line_range[0]), int(line_range[1]) + 1))

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
# Find intersects between two ranges (a and b)
###
def intersects(a, b):
    overlap_ranges = []
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

    for i in range(len(overlap_ranges) - 1):
        if overlap_ranges[i][1] in (overlap_ranges[i + 1][0], overlap_ranges[i + 1][0] - 1):
            overlap_ranges[i:i + 2] = [[overlap_ranges[i][0], overlap_ranges[i + 1][1]]]

    return overlap_ranges


###
# Extract the factors overlaps in percentage
# Inputs:
#   - genomes_directory : directory of the genomes of species
#   - main_proxies_directory : main directory of the proxies obtained through the factor_proxies script
#   - all_factors : list of all wanted factors
#   - species : name of the wanted species
# Output:
#   - A list which will contain all pairwise overlap between the factors
###
def extract_overlap_percentages(genomes_directory, main_proxies_directory, all_factors, species):
    # Fetch this species records for the lengths
    records = fetch_fasta(extract_path(genomes_directory, species + '*')[0])

    # Directory of the different factor proxies of this species
    species_proxies_directory = main_proxies_directory + '/' + species

    species_overlaps = list()

    for each_factor in range(len(all_factors)):
        for each_comparison in range(len(all_factors)):
            # Avoid comparing factor overlap to himself...
            if each_factor == each_comparison:
                continue
            # We only do the "lower diagonal" -> don't want to do twice the same operation
            elif each_factor < each_comparison:
                factor_directory = species_proxies_directory + '/' + all_factors[each_factor]
                comparison_directory = species_proxies_directory + '/' + all_factors[each_comparison]

                # Catch any missing factor
                try:
                    # Find the sum of all nucleotide which are either factor or comparison factor
                    pure_factor_length = sum([count_pure_record_length(records[each_record], factor_directory)
                                              for each_record in range(len(records))])
                    pure_comparison_length = sum([count_pure_record_length(records[each_record], comparison_directory)
                                                  for each_record in range(len(records))])

                    # Make sure neither are empty
                    if pure_factor_length and pure_comparison_length:
                        # This vector will contain the sum of all overlapping nucleotides between factor and comparison
                        length_overlapping = 0

                        for each_record in range(len(records)):
                            factor_ranges = extract_ranges(factor_directory + '/' + records[each_record].id)
                            comparison_ranges = extract_ranges(comparison_directory + '/' + records[each_record].id)
                            # Find all overlapping ranges
                            overlap_ranges = intersects(factor_ranges, comparison_ranges)

                            for each_range in overlap_ranges:
                                # Add the range length
                                length_overlapping += len(range(each_range[0], each_range[1] + 1))

                        factor_overlap = (length_overlapping / pure_factor_length) * 100
                        comparison_overlap = (length_overlapping / pure_comparison_length) * 100

                        # We thus have: the species/the factor we look at at/
                        # the % of this factor which overlap/the comparison factor
                        species_overlaps.append([species, all_factors[each_factor],
                                                 str(math.floor(factor_overlap * 100) / 100),
                                                 all_factors[each_comparison]])
                        species_overlaps.append([species, all_factors[each_comparison],
                                                 str(str(math.floor(comparison_overlap * 100) / 100)),
                                                 all_factors[each_factor]])
                # No factor in this species
                except:
                    # We need some dummy answers for the plots later on
                    species_overlaps.append([species, all_factors[each_factor],
                                             str(0), all_factors[each_comparison]])
                    species_overlaps.append([species, all_factors[each_comparison],
                                             str(0), all_factors[each_factor]])

    return species_overlaps


all_species = [os.path.basename(each) for each in extract_path("config/species/", '*')]
all_factors = [os.path.basename(each) for each in extract_path("config/factors/", '*')]

# Directory which lead to the different species proxies
main_proxies_directory = "../files/factor_proxies/" + str(window_size)

genomes_directory = "../data/genomes/"

# This list will keep all the factors percentages, in all species
all_percentages = Parallel(n_jobs=n_threads)(delayed(extract_overlap_percentages)
                                             (genomes_directory, main_proxies_directory, all_factors,
                                              all_species[each_species])
                                             for each_species in range(len(all_species)))

# The percentages may then be stored
with open(output, 'w') as outfile:
    for each_species in all_percentages:
        for each_perc in each_species:
            outfile.write('\t'.join(each_perc) + '\n')
