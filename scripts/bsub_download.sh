#!/usr/bin/env bash

# To use : bsub < ./bsub_download.sh
# Will use multithread (-n = number of cores, -R = same host)
# Takes about 3 minutes to complete


#BSUB -L /bin/bash
#BSUB -e error_dl.txt
#BSUB -J dl
#BSUB -n 4
#BSUB –R "span[ptile=4]"
#BSUB -M 10485760

bash download_genomes.sh
