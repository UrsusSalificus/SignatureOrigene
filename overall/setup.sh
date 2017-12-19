#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# Will setup the whole analysis through user commands
# Usage : bash setup.sh other_snakemake_arguments
# ASCII arty letters created through    http://patorjk.com/software/taag/#p=display&f=Chunky&t=Example

# Any argument can be passed to snakemake as argument of this setup script
# E.g. number of cores (--cores #) or dry run (-np = will just show the jobs to do, not really running them)
snakemake_arguments="$*"

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
echo "You are about to perform data health checks"
echo "________________________________________________________________________________"
echo ""
read -n 1 -s -r -p "Press any key to continue"
echo ""
echo "________________________________________________________________________________"
echo "

"

########### SOURCING THE CHOOSING FUNCTIONS ###########
source ../scripts/setup_functions.sh

########### FIGURES ###########
title="config/temp/title.txt"
cat << "EOF" > $title
 _______ __
|    ___|__|.-----.--.--.----.-----.-----.
|    ___|  ||  _  |  |  |   _|  -__|__ --|
|___|   |__||___  |_____|__| |_____|_____|
            |_____|

EOF
intro="Please chose among the following list which figures should be included in the analysis."
table="../input/overall_specific/table_figures.txt"
choice=fixed
good_inputs="../input/overall_specific/good_figures.txt"
abbrev="../input/overall_specific/abbrev_figures.txt"
out_dir="config/figures/"
nice="../input/overall_specific/all_figures.txt"

confirm



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



########### FACTORS ###########
title="config/temp/title.txt"
cat << "EOF" > $title
 _______              __
|    ___|.---.-.----.|  |_.-----.----.-----.
|    ___||  _  |  __||   _|  _  |   _|__ --|
|___|    |___._|____||____|_____|__| |_____|

EOF
intro="Please chose among the following list which factor should be included in the analysis."
table="../input/purifying_specific/table_factors.txt"
choice=fixed
good_inputs="../input/purifying_specific/good_factors.txt"
abbrev="../input/purifying_specific/abbrev_factors.txt"
out_dir="config/factors/"
nice="../input/purifying_specific/all_factors.txt"

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
cd ../factors
FACTORS=$( find * )
cd ../figures
FIGURES=$( find * )
cd ../..

# Remember where we are :
go_back=$( pwd )

# The first part is to extract and clean factor ranges:
for each_species in $SPECIES; do
    for each_factor in $FACTORS; do
        # If feature = UTR, do not compute for S. cerevisiae and E. coli
        if [[ $each_factor == 'UTR' ]] && \
            ([[ $each_species == 's_cerevisiae' ]] || [[ $each_species == 'e_coli' ]]); then
            echo "$each_factor not computed for $each_species"
        elif [[ $each_factor == 'intron' && $each_species == 'e_coli' ]] ; then
            echo "$each_factor not computed for $each_species"
        else
            # Finding factor ranges
            snakemake $snakemake_arguments \
                ../data/following/factor_proxies/$each_species/$each_factor\_proxies_done.txt
        fi
    done
    for each_window in $WINDOWS; do
        # After finding all the factors' ranges, we must clean them from overlaps
        snakemake $snakemake_arguments \
            ../data/following/factor_filtered/$each_window/$each_species\_done.txt
    done
done

# Second part will extract and plot the various percentages
for each_window in $WINDOWS; do
    for each_figure in $FIGURES; do
        snakemake $snakemake_arguments files/results/$each_window/$each_figure\_all_species.png
    done
done