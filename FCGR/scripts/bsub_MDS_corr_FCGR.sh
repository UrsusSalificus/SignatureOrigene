#!/usr/bin/env bash

# Usage (in the FCGR directory) : bsub < scripts/bsub_MDS_corr_FCGR.sh
# Will use multithread (-n = number of cores, -R = same host)
# Will use 15 Gb, on the same host

#BSUB -L /bin/bash
#BSUB -e error_MDS_cor_F.txt
#BSUB -J MDS_cor_F
#BSUB -n 16
#BSUB –R "span[hosts=1]"
#BSUB -M 256000000
#BSUB –R "rusage[mem=250000]"

module add R/3.3.2
module add Utility/snakemake/3.11.2;

bash cluster_setup.sh --cores 14