"""This script will clean the genome fasta files to only have wanted records
(only full chromosomes were kept for simplicity)
"""
from Bio import SeqIO
import sys

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Species genome path:
species_genome = str(sys.argv[1])
# Species abbreviation:
species = '_'.join(str(species_genome.split('/')[-1]).split('_')[:2])


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

# Store all the records we which to keep
with open("input/wanted_records.txt") as records_ids:
    do_want = [each_id.split()[0].strip() for each_id in records_ids]

records = fetch_fasta(species_genome)
cleaned_records = [each_record for each_record in records if each_record.id in do_want]

with open(species_genome, "w") as cleaned_output:
    SeqIO.write(cleaned_records, cleaned_output, "fasta")