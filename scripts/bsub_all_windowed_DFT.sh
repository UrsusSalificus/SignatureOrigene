#!/usr/bin/env bash

# To use : bsub < ./bsub_all_windowed_DFT.sh
# Will use multithread (-n = number of cores, -R = same host)
# Will use 15 Gb, on the same host

#BSUB -L /bin/bash
#BSUB -e error_DFT_a.txt
#BSUB -J DFT_a
#BSUB -n 32
#BSUB –R "span[ptile=32]"
#BSUB -M 256000000
#BSUB –R "rusage[mem=250000]"

python3 all_windowed_DFT.py 150000 16