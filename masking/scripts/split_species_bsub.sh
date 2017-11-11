#!/usr/bin/env bash

SPECIES="hsap_sample mmus_sample c_elegans d_melanogaster a_thaliana s_cerevisiae e_coli"
for each_species in $SPECIES; do
  echo """
  #!/usr/bin/env bash

  # Usage (in the masking directory) : bsub < scripts/bsub_masking.sh
  # Will use multithread (-n = number of cores, -R = same host)
  # Will use 16 Gb, on the same host

  #BSUB -L /bin/bash
  #BSUB -e error_masking_$each_species.txt
  #BSUB -J masking_$each_species
  #BSUB -n 8
  #BSUB -R "span[hosts=1]"
  #BSUB -M 16384000
  #BSUB -R "rusage[mem=16000]"

  module add R/3.3.2
  module add Utility/snakemake/3.11.2;

  bash cluster_setup.sh $each_species --cores 6
  """ > temp_bsub.sh

  bsub < temp_bsub.sh

  rm temp_bsub.sh
done
