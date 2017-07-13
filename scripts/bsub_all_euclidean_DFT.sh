#!/usr/bin/env bash
# To use : bsub < ./bsub_all_euclidean_DFT.sh
# Will use multithread (-n = number of cores, -R = same host)
# Will use 15 Gb, on the same host

#BSUB -L /bin/bash
#BSUB -e error_euc_a.txt
#BSUB -J euc_a
#BSUB -n 8
#BSUB –R "span[ptile=8]"
#BSUB -M 15728640
#BSUB –R "rusage[mem=15360]"

python3 all_Euclidean.py 150000 3 DFTs