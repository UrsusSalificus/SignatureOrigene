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
echo "You are about to mask the genome to only have certain features. Both these new"
echo "sequences (without the feature) and pure factor sequences (with only the factor)"
echo "will then be compared to the representative of the whole genome's kmer"
echo "frequencies (called the center)."
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
Type sample size (eg. 500 1000 for both 500 and 1000 windows per species), then press
[ENTER]:
Note: if a species genome is not big enough for the wanted number of sample windows
the whole genome is taken instead.
EOF
choice=free
good_inputs=1234567890
abbrev="config/temp/abbrev.txt"
cat << "EOF" > $abbrev
free
EOF
out_dir="config/samples/"
nice="config/temp/nice.txt"
cat << "EOF" > $nice
Windows per species =
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
cd ../samples
SAMPLES=$( find * )
cd ../windows
WINDOWS=$( find * )
cd ../factors
FACTORS=$( find * )
cd ../kmer
KMER=$( find * )
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
    # After finding all the factors' ranges, we must clean them from overlaps
    snakemake $snakemake_arguments \
        ../data/following/factor_proxies/uncategorized/$each_species\_done.txt
done


# The second part goes from these ranges to find sequences where the factor is absent
for each_species in $SPECIES; do
    for each_sample in $SAMPLES; do
        for each_window in $WINDOWS; do
            for each_kmer in $KMER; do
                # A) We need the whole genome distance matrix, which will be computed through the scaling snakemake
                cd ../scaling
                snakemake $snakemake_arguments \
                    files/distances/manhattan/$each_window\_$each_sample\_$each_kmer/$each_species\_dist_matrix.RData
                cd $go_back

                # B) We now will compute both the masked and pure factors distance matrix VS center
                for each_factor in $FACTORS; do
                    # If feature = UTR, do not compute for S. cerevisiae and E. coli
                    if [[ $each_factor == 'UTR' ]] && \
                        ([[ $each_species == 's_cerevisiae' ]] || [[ $each_species == 'e_coli' ]]); then
                        echo "$each_factor not computed for $each_species"
                    elif [[ $each_factor == 'intron' && $each_species == 'e_coli' ]] ; then
                        echo "$each_factor not computed for $each_species"
                    else

                        # B.1) Masked will be done all on the masking directory
                        snakemake $snakemake_arguments \
                            files/distances/manhattan/$each_window\_$each_sample\_$each_kmer/$each_species\_$each_factor\_masked_vs_center_dist_matrix.RData

                        # B.2) Pure will need the FCGRs found in the purifying directory
                        cd ../purifying
                        snakemake $snakemake_arguments \
                            files/FCGRs/$each_window\_$each_sample\_$each_kmer/$each_species\_$each_factor\_pure_FCGRs.txt
                        cd $go_back
                        # B.3) Which will then be used to find the pure VS center distance matrix
                        snakemake $snakemake_arguments \
                            files/distances/manhattan/$each_window\_$each_sample\_$each_kmer/$each_species\_$each_factor\_pure_vs_center_dist_matrix.RData
                    fi
                done

                # C) Our control will be the whole genome VS center distance
                snakemake $snakemake_arguments \
                        files/distances/manhattan/$each_window\_$each_sample\_$each_kmer/$each_species\_whole_vs_center_dist_matrix.RData
                # TODO: find why we have to recompute the distance this way instead of using the whole distance matrix

                # D) We would finally be able to use both distance matrices to produces boxplots
                snakemake $snakemake_arguments \
                        files/results/$each_window\_$each_sample\_$each_kmer/$each_species\_boxplots_all_factors.png
            done
        done
    done
done
