#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# Will setup the random analysis through user commands
# Usage : bash setup.sh other_snakemake_arguments
# ASCII arty letters created through    http://patorjk.com/software/taag/#p=display&f=Chunky&t=Example

# The number of cores used can be set as the first argument of the script
n_cores=$1
while [[ "$n_cores" == "" ]]; do
  echo "Please, enter the number of cores you wish to use (e.g. 2)."
  read input
  n_cores=$input
done

mkdir -p config/temp

echo " _______ __                     __
|     __|__|.-----.-----.---.-.|  |_.--.--.----.-----.
|__     |  ||  _  |     |  _  ||   _|  |  |   _|  -__|
|_______|__||___  |__|__|___._||____|_____|__| |_____|
            |_____|                                   "
echo " _______        __
|       |.----.|__|.-----.-----.-----.-----.
|   -   ||   _||  ||  _  |  -__|     |  -__|
|_______||__|  |__||___  |_____|__|__|_____|
                   |_____|                  "
echo "________________________________________________________________________________"
echo ""
echo "Welcome to the setup of the analysis!"
echo "You are about to perform a comparison of the genomic signatures of n random     "
echo "sequences per species."
echo "Note: computational time does not matter anymore when taking species with big "
echo "genomes (as long as n is < 100)."
echo "________________________________________________________________________________"
echo ""
read -n 1 -s -r -p "Press any key to continue"
echo ""
echo "________________________________________________________________________________"
echo "

"

########### SOURCING THE CHOOSING FUNCTIONS ###########
source ../scripts/setup_functions.sh

########### SPECIES ###########
title="config/temp/title.txt"
cat << "EOF" > $title
     _______                    __
    |     __|.-----.-----.----.|__|.-----.-----.
    |__     ||  _  |  -__|  __||  ||  -__|__ --|
    |_______||   __|_____|____||__||_____|_____|
             |__|

EOF
intro="Please chose among the following list which species should be included in the analysis."
table="../input/table_species.txt"
choice=fixed
good_inputs="../input/good_species.txt"
abbrev="../input/abbrev_species.txt"
out_dir="config/species/"
nice="../input/all_species.txt"

confirm

# Check there is at least two species...
if [[ $( ls config/species | wc -w ) -lt 2 ]] ; then
    echo "Cannot compute the analysis with only one species, please select more!"
    exit
fi

########### WINDOW SIZE TYPES ###########
title="config/temp/title.txt"
cat << "EOF" > $title
     ________ __           __                             __
    |  |  |  |__|.-----.--|  |.-----.--.--.--.    .-----.|__|.-----.-----.
    |  |  |  |  ||     |  _  ||  _  |  |  |  |    |__ --||  ||-- __|  -__|
    |________|__||__|__|_____||_____|________|    |_____||__||_____|_____|

EOF
intro="Please chose the size of the windows on which will be computed the Chaos Game \nRepresentation"
table="config/temp/table.txt"
cat << "EOF" > $table
Type the window size (in base pairs, eg. 5000 150000 for both 5kbp and 15kb
windows), then press [ENTER]:
Here are displayed the minimum window size, for each to have at least a chance
to see each k-mer once:

    Window size | k-mer |
        192     |   3   |
        1024    |   4   |
        5120    |   5   |
        24576   |   6   |
        114688  |   7   |
        524288  |   8   |
        2359296 |   9   |

5kb and 15kb were initially chosen for the analysis, with k-mer of respectively
k=4 and k=7
NOTE: if you were to chose FCGR as genomic signature's type later on, each set
of window size/k-mer size would be computed, heavily increasing the
computational time. It this advised to chose one set per analysis.
EOF
choice=free
good_inputs=1234567890
abbrev="config/temp/abbrev.txt"
cat << "EOF" > $abbrev
free
EOF
out_dir="config/windows/"
nice="config/temp/nice.txt"
cat << "EOF" > $nice
Windows of size (in base pairs) :
EOF

confirm



########### SAMPLE SIZE ###########
title="config/temp/title.txt"
cat << "EOF" > $title
     _______                        __                    __
    |     __|.---.-.--------.-----.|  |.-----.    .-----.|__|.-----.-----.
    |__     ||  _  |        |  _  ||  ||  -__|    |__ --||  ||-- __|  -__|
    |_______||___._|__|__|__|   __||__||_____|    |_____||__||_____|_____|
                            |__|
EOF
intro="Please chose the sample size (number of windows per species)"
table="config/temp/table.txt"
cat << "EOF" > $table
________________________________________________________________________________
Type sample size (eg. 10 20 for both 10 and 20 windows per species), then press
[ENTER]:
EOF
choice=free
good_inputs=1234567890
abbrev="config/temp/abbrev.txt"
cat << "EOF" > $abbrev
free
EOF
out_dir="config/sample/"
nice="config/temp/nice.txt"
cat << "EOF" > $nice
Windows per species =
EOF

confirm



########### K-MER SIZE ###########
title="config/temp/title.txt"
cat << "EOF" > $title
     __  __                                       __
    |  |/  |_____.--------.-----.----.    .-----.|__|.-----.-----.
    |     <______|        |  -__|   _|    |__ --||  ||-- __|  -__|
    |__|\__|     |__|__|__|_____|__|      |_____||__||_____|_____|

EOF
intro="Please chose the size of the word (k-mer) for the FCGR computation"
table="config/temp/table.txt"
cat << "EOF" > $table
________________________________________________________________________________
Type k size (eg. 4 7 for both k=4 and k=7), then press [ENTER]:
Here are displayed the minimum window size, for each to have at least a chance
to see each k-mer once:

    Window size | k-mer |
        4       |   1   |
        32      |   2   |
        192     |   3   |
        1024    |   4   |
        5120    |   5   |
        24576   |   6   |
        114688  |   7   |
        524288  |   8   |
        2359296 |   9   |

5kb and 15kb were initially chosen for the analysis, with k-mer of respectively
k=4 and k=7
EOF
choice=free
good_inputs=1234567890
abbrev="config/temp/abbrev.txt"
cat << "EOF" > $abbrev
free
EOF
out_dir="config/kmer/"
nice="config/temp/nice.txt"
cat << "EOF" > $nice
k =
EOF

confirm




echo " ___ ___ _______ _______      ______ _______ ______ __  __ __
|   |   |       |   |   |    |   __ \       |      |  |/  |  |
 \     /|   -   |   |   |    |      <   -   |   ---|     <|__|
  |___| |_______|_______|    |___|__|_______|______|__|\__|__|
                                                              "

# Cleaning
cd config
rm -r temp
cd ..


# Do the entire analysis for each combination
cd config/species
SPECIES=$( find * )
cd ../windows
WINDOWS=$( find * )
cd ../kmer
KMER=$( find * )
cd ../sample
SAMPLE=$( find * )
cd ../..


for each_species in $SPECIES; do
    # We will have to check the downloaded files
    go_back=$( pwd )
    cd ..
    genome_file=data/genomes/$each_species\_genomes.fna
    feature_file=data/factors/features/$each_species\_feature_table.txt
    repeat_file=data/factors/repeats/$each_species\_repeats.txt
    # If any of those is missing, download again
    if [[ ! -f $genome_file || ! -f $feature_file || ! -f $repeat_file ]]; then
        bash scripts/download_genomes.sh $each_species $genome_file $feature_file $repeat_file
        # Clean of unwanted records
        python3 scripts/keep_wanted_records.py $genome_file
    fi
    cd $go_back
done

for each_sample in $SAMPLE; do
    for each_window in $WINDOWS; do
        for each_kmer in $KMER; do
            python3 scripts/random_GS.py $each_sample $each_window $each_kmer $n_cores
            Rscript scripts/random_MDS.R
        done
    done
done

# Clean everything
rm -r config
rm -r temp