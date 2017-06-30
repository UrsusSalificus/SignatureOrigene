#!/usr/bin/env bash

# To use : bsub < ./bsub_CGR_whole.sh
# Will use multithread on 2 hosts (-n = number of cores, -R = same host)
# Takes about 3 minutes to complete


#BSUB -L /bin/bash
#BSUB -e error_CGR_w.txt
#BSUB -J CGR_w
#BSUB -n 4
#BSUB –R "span[ptile=4]"
#BSUB -M 10485760

bash CGR_whole_all.sh