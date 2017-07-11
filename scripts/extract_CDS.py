# This script concatenate all the different FCGRs in a single file
import CGR_functions as fn
import math


# Wanted window size:
#window_size = int(sys.argv[1])
window_size = 150000
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
    species_CDS = '/'.join(['../files/features/CDS', window_in_kb, species[each_species] + '_CDS'])
    fn.checking_parent(species_CDS)

    # Fetch all the records from this species fasta
    records = fn.fetch_fasta(species_genome)

    with open(species_table, 'r') as feature_table, open(species_CDS, 'w') as outfile:
        feature_table.readline()
        actual_line = feature_table.readline().split('\t')
        for each_record in range(len(records)):
            # Each element of this list represents a nucleotide
            record_proxy = [0] * len(records[each_record].seq)
            while records[each_record].id == actual_line[6]:
                # If it is a CDS
                if actual_line[0] == 'CDS':
                    # For each nucleotide from start to end of the CDS:
                    for each_nucleotide in range(int(actual_line[7]), int(actual_line[8])):
                        record_proxy[each_nucleotide] = 1
                actual_line = feature_table.readline().split('\t')
                # If we get at the last line, actual_line only have one empty entry, which can be detected by
                # calling the second element
                try:
                    actual_line[1]
                except:
                    break
            if len(records[each_record].seq) > window_size:
                n_windows = math.floor(len(records[each_record].seq) / window_size)  # Number of windows
                for start in range(0, n_windows):
                    window = str(records[each_record].seq)[(start * window_size):((start * window_size) + window_size)]
                    # If any character in the sequence is NOT a standard nucleotides (including unknown nucleotides),
                    # do NOT compute:
                    if not any([c not in 'ATCGatcg' for c in window]):
                        region_CDS = sum(record_proxy[(start * window_size):((start * window_size) + window_size)])
                        outfile.write(records[each_record].id + '\t')
                        outfile.write(str(region_CDS/window_size) + '\n')
