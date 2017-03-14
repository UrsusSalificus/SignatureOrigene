import itertools as it
import subprocess as sb
import time
import random
from Bio import SeqIO

class Sequence:
    def __init__(self, tetramers={}, records):
        self.tetramers = tetramers
        self.records = records

    def read_fasta(self, fasta_file):
        # Will only take the first sequence of the fasta file
        try:
            self.records = list(SeqIO.parse(fasta_file, "fasta"))[0]
        except:
            print("Cannot open %s" % fasta_file)
            quit()
        if len(list(SeqIO.parse(fasta_file, "fasta"))) > 1:
            print("Multiple sequences, will only compute CGR for the first one")

    def tetramer_frequency_whole(self, window_size=5000):
        # Compute each possible tetramer using it.product
        for i in it.product('ATCG', repeat=4):
            # For each tetramer, input an entry in dictionary, with 0 as sister value.
            # Use ''.join to concatenate the it.product object from a list of individual strings to single string
            self.tetramers[''.join(i)] = 0

        # Will now look for each of these tetramer
        for product in self.tetramers:
            self.tetramers[product] = self.records.seq.count(product) / (len(self.records.seq)-3)

    def tetramer_frequency_moving_window(self, window_size=5000):
        max_n_windows=0.75*len(self.seq)/window_size
        whole_freq=

        for window_number in range(0, max_n_windows):
            start=random.randrange(0, len(self.seq)-5000)
            # TODO: finish this !
            end=start+5000

            # Compute each possible tetramer using it.product
            for i in it.product('ATCG', repeat=4):
                # For each tetramer, input an entry in dictionary, with 0 as sister value.
                # Use ''.join to concatenate the it.product object from a list of individual strings to single string
                self.tetramers[''.join(i)] = 0

            for product in self.tetramers:
                self.tetramers[product] = self.seq[start:end].count(product) / window_size



for seq in SeqIO.parse("opsin_Amphiprion.fasta", "fasta"):
    print(len(seq.seq))
    sequences=seq.seq
print(sequences)

records = SeqIO.parse("opsin_Amphiprion.fasta", "fasta")
for each in records:
    print each.id
sizes = len([len(rec) for rec in records])
records[0].seq.count("AT")


print(random.randrange(0, 10000-5000))
0.75*100000/5000
for lol in range(0, 2):
    print(lol)

f = open('test')
seq = f.readline()
seq[0:2]

timer_start = time.time()

test_seq = Sequence()
test_seq.read_from_concatanated('dmel_all_chromosome_single_seq')

timer_end = time.time() - timer_start

print( timer_end)

timer_start = time.time()
test_seq.tetramer_frequency()
timer_end = time.time() - timer_start
print( timer_end)

test_seq.tetramers

