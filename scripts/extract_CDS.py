# This script concatenate all the different FCGRs in a single file
import CGR_functions as fn
from Bio import SeqIO
import glob
import math
import sys
import os

# Wanted window size:
window_size = int(sys.argv[1])
window_in_kb = str(window_size)[:-3] + 'kb'

# Get all the species abbreviation for the study
species = fn.get_species()

CDS_all_records = []
for each_species in range(len(species)):
    # Due to non-consistent pattern of file name, whole genome ('genomic') is used in multiple
    # fasta files (CDS or RNA only, and any nucleotides) names. We must thus reconstruct the exact path.
    all_files = fn.extract_path('../data/genomes/', str(species[each_species] + '*'))
    index = all_files[0].split('/')[3].split('_')[2] + '_' + all_files[0].split('/')[3].split('_')[3]
    pattern_genome = str(species[each_species] + '_' + index + '_genomic*')
    species_genome = fn.extract_path('../data/genomes/', pattern_genome)[0]

    # For the feature table, it is much easier (only one per species)
    species_table = fn.extract_path('../data/genomes/' + species[each_species], '*_feature_table*')[0]

    # Finally, the output file :
    # TODO : you have to find a way to stock everything
    species_CDS = '/'.join(['../files/CDS/', window_in_kb])

    records = fn.fetch_fasta(species_genome)
    for each_record in range(len(records)):
        # Each element of this list represents a nucleotide
        record_proxy = [0] * len(records[each_record].seq)
        if len(records[each_record].seq) > window_size:
            n_windows = math.floor(len(records[each_record].seq) / window_size)  # Number of windows
            windows_kept = []
            for start in range(0, n_windows):
                window = str(records[each_record].seq)[(start * window_size):((start * window_size) + window_size)]
                # If any character in the sequence is NOT a standard nucleotides (including unknown nucleotides), do NOT compute:
                if not any([c not in 'ATCGatcg' for c in window]):
                    windows_kept.append(start)
                    with open(species_table, 'r') as feature_table:
                        for each_line in feature_table:
                            split_line = each_line.split('\t')
                            if split_line[6] == records[each_record].id and split_line[0] == 'CDS':
                                for each_nucleotide in range(int(split_line[7]), int(split_line[8])):
                                    record_proxy[each_nucleotide] = 1


for each_species in range(len(species)):
    # For the feature table, it is much easier (only one per species)
    species_table = fn.extract_path('../data/genomes/' + species[each_species], '*_feature_table*')[0]
    with open(species_table, 'r') as feature_table:
        feature_table.readline()
        print(feature_table.readline().split('\t')[7])
        print(feature_table.readline().split('\t')[8])





    with open(species_table, 'r') as feature_table:
        for each_record in range(len(all_records)):

            for each_line in feature_table:
                split_line = each_line.split()
                if split_line[5] == :
                    if split_line[0] == 'CDS':
                        starts_I.append(split_line[7])
                        ends_I.append(split_line[8])


            all_region = fn.extract_path(all_records[each_record] + '/', '*')





            all_regions = fn.extract_path(each_record+ '/', '*')
            record_name = each_record.split('/')[-1]
            for each_region in range(len(all_regions)):
                outfile.write

                CGR_directory = '/'.join(['../files/CGRs', window_in_kb, species[each_species]]) + '/'
                all_records = fn.extract_path(CGR_directory + '/', '*')
                records_names = [each_record.split('/')[-1] for each_record in all_records]