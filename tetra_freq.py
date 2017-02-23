import itertools as it

# Create empty dictionary that will contain all the 256 different possible tetramers
tetramer = {}

# Compute each possible tetramer using it.product
for i in it.product('ATCG', repeat=4):
    # For each tetramer, input an entry in dictionary, with 0 as sister value.
    # Use ''.join to concatenate the it.product object from a list of individual strings to single string
    tetramer[''.join(i)] = 0

tetramer

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


test_seq = Sequence()

test_seq.read_from_fasta('opsin_Amphiprion.fasta')
test_seq.tetramer_frequency()

test_seq.tetramers

