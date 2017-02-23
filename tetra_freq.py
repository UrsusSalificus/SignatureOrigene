import itertools as it

class Sequence:
    def __init__(self, seq='', tetramers={}):
        self.seq = seq
        self.tetramers = tetramers

        # Compute each possible tetramer using it.product
        for i in it.product('ATCG', repeat=4):
            # For each tetramer, input an entry in dictionary, with 0 as sister value.
            # Use ''.join to concatenate the it.product object from a list of individual strings to single string
            self.tetramers[''.join(i)] = 0

    def read_from_fasta(self, fasta_file):
        try:
            self.fasta = open(fasta_file, 'r')
        except:
            print("Cannot open %s" % fasta_file)
            quit()

        # Now, take only lines without > to have only a big sequence
        for lines in self.fasta:
            line = lines.rstrip()
            if '>' not in line:
                self.seq += line

    def tetramer_frequency(self):
        for product in self.tetramers:
            self.tetramers[product] = self.seq.count(product) / (len(self.seq)-3)


import time

timer_start = time.time()

test_seq = Sequence()
test_seq.read_from_fasta('dmel-2L-chromosome-r5.49.fasta')

timer_end = time.time() - timer_start

print( timer_end)


with open('dmel_2L_1sequence', 'w') as f:
    f.write(test_seq.seq)

test_seq.seq

timer_start = time.time()
test_seq.tetramer_frequency()
timer_end = time.time() - timer_start
print( timer_end)

test_seq.tetramers

