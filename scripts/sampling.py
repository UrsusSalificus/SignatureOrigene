#!/usr/bin/env python3

""" Extract the sample chromosome out of H. sapiens or M. musculus genome
"""
from Bio import SeqIO
import sys


__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Species genome path:
species_genome = str(sys.argv[1])

# Here we have both H. sapiens or M. musculus (respectively) chromosome 10 ids
chrom10_ids = ('NC_000010.11', 'NC_000076.6')

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

records = fetch_fasta(species_genome)

chrom10 = [each_record for each_record in records if each_record.id in chrom10_ids]

with open(species_genome, "w") as outfile:
    SeqIO.write(chrom10, outfile, "fasta")


all_chr10 = list()
with open('/home/titouan/PycharmProjects/Master/Main/FCGR/files/FCGRs/15000_7/h_sapiens_FCGRs.txt', 'r') as infile:
    for each_line in infile:
        line = each_line.split()
        if [line[0] == 'NC_000010.11']:
            all_chr10.append(line)


