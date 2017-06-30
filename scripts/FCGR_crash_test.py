# Making sure there is no mistakes in the FCGR_from_CGR function
import CGR_functions
import itertools as it
from Bio import SeqIO

# This script enable to make sure the nucleotide frequencies obtained through own FCGR function
# AND Python's one by one dictionary checking function are identical
# Beta globin cds of human is used to perform this test
fasta_file = '../data/h_sapiens_globin.fasta'

records = list(SeqIO.parse(fasta_file, "fasta"))[0]

# Chose k-mer size:
k_size = 2


# Reference: checking every position, one at a time:
def one_by_one(records, k_size):
    k_mer = {}
    for i in it.product('ATCG', repeat=k_size):
        # For each tetramer, input an entry in dictionary, with 0 as sister value.
        # Use ''.join to concatenate the it.product object from a list of individual strings to single string
        k_mer[''.join(i)] = 0
    for product in k_mer:
        for each_base in range(len(records) - k_size + 1):
            # since some genomes contain cap sequences 'N' check whether key exists and add 1 if it does
            window = str(records.seq)[each_base:each_base + k_size]
            if window == product:
                k_mer[window] += 1
    return(k_mer)
ref = one_by_one(records, k_size)

# The FCGR function:
CGR = CGR_functions.CGR_coordinates(str(records.seq), '')   # '' -> don't want an outfile
FCGR = CGR_functions.FCGR_from_CGR(k_size, CGR, '')   # '' -> don't want an outfile
FCGR_labels = CGR_functions.FCGR_indexes(k_size)


# Now compare:
if sum(FCGR) == sum(ref.values()):
    print('Same length!')

errors = []
for number_k_mer in range(len(FCGR_labels)):
    each_k_mer = FCGR_labels[number_k_mer]
    if not FCGR[number_k_mer] == ref[each_k_mer]:
        errors.append(str('Problem with ' + each_k_mer + ': count = ' + str(k_mer[each_k_mer]) +
              ' / FCGR = ' + str(FCGR[number_k_mer])))
if not errors:
    print('No error detected!')
else:
    for each_error in errors:
        print(each_error)
