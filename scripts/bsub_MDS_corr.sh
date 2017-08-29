#!/usr/bin/env bash

# To use : bsub < ./bsub_all_windowed_DFT.sh
# Will use multithread (-n = number of cores, -R = same host)
# Will use 15 Gb, on the same host

#BSUB -L /bin/bash
#BSUB -e error_MDS_cor.txt
#BSUB -J MDS_cor
#BSUB -n 32
#BSUB –R "span[hosts=2]"
#BSUB -M 256000000
#BSUB –R "rusage[mem=250000]"

bash ./DFT/cluster_setup.sh --cores 14
bash ./FCGR/cluster_setup.sh --cores 14