#!/usr/bin/env bash

# Designed only for concatenation of fasta containing a single sequence (only one >).

# How to use this :
# 1) Launch bash terminal
# 2) cd (change directory) to the directory/folder containing your fasta file you want to concatenate
# (optional) Link to this directory this bash file using the following command line:
#               ln -s path_to_concatenate_fasta.sh
#       - Where path is your path to the script file (ex. ./Scripts/concatenate_fasta.sh)
# 3) Run the script using the following command line :
#               sh path_to_concatenate_fasta.sh fasta_file.fasta concatenated_name
#       - Where your replace path_to_concatenate_fasta.sh by the path to this script file
#       (eg. ./Scripts/concatenate_fasta.sh ; or just write concatenate_fasta.sh if did optional step),
#       fasta_file.fasta by your fasta file name, and concatenated_name by your desired output file name.

### Delete the line containing > ; then concatenate by deleting \n
# sed : stream editor
# ^ : will anchor the following pattern to the start
# d : delete
sed '/^>/d' $1 | tr -d '\n' > $2


