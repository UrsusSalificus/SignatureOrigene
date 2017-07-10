#!/usr/bin/env bash

# To use : bsub < ./bsub_all_windowed_FCGR.sh
# Will use multithread (-n = number of cores, -R = same host)
# Will use 15 Gb, on the same host

#BSUB -L /bin/bash
#BSUB -e error_FCGR_a.txt
#BSUB -J FCGR_a
#BSUB -n 32
#BSUB –R "span[ptile=32]"
#BSUB -M 102400000
#BSUB –R "rusage[mem=100000]"

python3 all_windowed_FCGR.py 9 150000 16