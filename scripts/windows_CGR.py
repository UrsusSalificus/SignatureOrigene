# This script selects the function required for computing windowed FCGRs of a sequence.
import CGR_functions
from joblib import Parallel, delayed
import math
import subprocess

# Path to the file containing the sequence one wants the CGR computed on
fasta_file = '../data/genomes/e_coli/GCF_000023665.1_ASM2366v1_genomic.fna'
# Path to the output file, which will contain the x/y CGR coordinates
# If empty, will print the coordinates instead of writing a file.
CGR_outfile = '../data/genomes/e_coli/all_CGRs/CGR_region_'
# Path to the file which will be used in R to compute PCA
FCGRs_for_R = '../files/FCGRs/FCGRs_e_coli'
# Path to the output picture of PC1/PC2
PCA = '../files/PCA/PCA_e_coli'
# Size of the sliding window (not overlapping)
window_size = 150000
window_in_kb = '15kb'
# K-mer size
k_size = 3


def windowed_FCGR(fasta_file, CGR_outfile, FCGRs_for_R, window_size, window_in_kb, k_size, PCA):
    # First step: checking that all the parents directory are here:
    all_paths = [CGR_outfile, FCGRs_for_R, PCA]
    for each_path in all_paths:
        CGR_functions.checking_parent(each_path)


    # After fetching the cleaned sequence, we will see how many windows we can fit (n_windows)
    cleaned_seq = CGR_functions.fetch_and_clean_fasta(fasta_file)
    n_windows = math.floor(len(cleaned_seq) / window_size)
    rounded_seq_len = window_size * n_windows

    # For each window, we want to compute the CGR
    # We use joblib.Parallel to reduce computation time:
    # n_jobs = number of threads -> here put the maximum number of cores that should be working on this
    # delayed will make tupple oof the form (function) (arguments of the function)
    Parallel(n_jobs=3)(delayed(CGR_functions.CGR_coordinates)
                       (str(cleaned_seq)[(start * window_size):((start * window_size) + window_size)],
                        (CGR_outfile + str(start)))
                       for start in range(0, n_windows))

    # We will now compute the k-mer Frequencies of the Chaos Game Representation (FCGR) of tall these CGR
    # Path to the CGr files:
    CGR_files = []
    for each_region in range(0, n_windows):
        CGR_files.append(CGR_outfile + str(each_region))

    # Will stock all the FCGRs in a list
    FCGRs = Parallel(n_jobs=3)(delayed(CGR_functions.FCGR_from_CGR)
                               (k_size, CGR_files[each_region], '')  # '' to have the ouput of the function as a list
                               for each_region in range(0, n_windows))

    indexes_of_FCGR = CGR_functions.FCGR_indexes(k_size)

    # Will stock FCGR as file easily readable for R
    with open(FCGRs_for_R, 'w') as outfile:
        for each_kmer in range(len(indexes_of_FCGR)):
            if not each_kmer == len(indexes_of_FCGR) - 1:
                outfile.write(indexes_of_FCGR[each_kmer] + '\t')
            else:
                outfile.write(indexes_of_FCGR[each_kmer])
        outfile.write('\n')
        for each_region in range(len(FCGRs)):
            for each_kmer in range(len(indexes_of_FCGR)):
                if not each_kmer == len(indexes_of_FCGR) - 1:
                    outfile.write(str(FCGRs[each_region][each_kmer]) + '\t')
                else:
                    outfile.write(str(FCGRs[each_region][each_kmer]))
            outfile.write('\n')

    # Now compute the PCA analysis on R
    subprocess.call(['Rscript', 'FCGR_PCA.R', FCGRs_for_R, window_in_kb, PCA])

    # Done ! Picture of PC1/PC2 available on @PCA path


windowed_FCGR(fasta_file, CGR_outfile, FCGRs_for_R, window_size, window_in_kb, k_size, PCA)
