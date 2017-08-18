#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# Will setup the whole analysis through user commands
# ASCII arty letters created through    http://patorjk.com/software/taag/#p=display&f=Chunky&t=Example

if [ -d config/* ]; then
    rm -r config/*
fi
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
echo "________________________________________________________________________________"
echo ""
#read -n 1 -s -r -p "Press any key to continue"
echo ""
echo "________________________________________________________________________________"
echo "

"

########### DEFINITION OF THE CHOOSING FUNCTION ###########
setup () {
    cat $1
    echo -e $2
    read -n 1 -s -r -p "Press any key to continue"
    echo ""
    echo "                  ### To do so ###"
    cat $3
    read input

    # Checking not empty
    while [[ "$input" == "" ]]; do
      echo "Please, enter at least one input."
      read input
    done

    if [[ "$4" == "fixed" ]]; then
        split_input=$( echo $input | grep -o . )
    fi

    # Checking this is the right input
    good_to_go=0
    good_input=$( echo $5 | grep -o . )
    until [[ "$good_to_go" -eq "1" ]]; do
        good_to_go=1
        for each_input in $split_input ; do
            # If not a letter from the list, stop the for loop, and must retype
            if [[ ! $good_input =~ (^|[[:space:]])"$each_input"($|[[:space:]]) ]]; then
                good_to_go=0
                continue
            fi
        done
        if [[ "$good_to_go" -eq "0" ]]; then
            echo "Please, enter a valid list (only numbers, in the right range, and at least one)"
            read input
            split_input=$( echo $input | grep -o . )
        fi
    done

    # If everything is alright, create the empty files with the right names
    # As we windows may be set freely, no hash tables are required,
    # we thus differentiate when there is fixed abbreviation, or when it is free choice.
    if [[ "$4" == "fixed" ]]; then
        declare -a all_input=($( cat $6 ))
        for each_input in $split_input; do
            if [ ! -f $7${all_input[$each_input]} ]; then
                mkdir -p $7 && touch $7${all_input[$each_input]}
            fi
        done
    else
        for each_input in $input; do
            if [ ! -f $7$each_input ]; then
                mkdir -p $7 && touch $7$each_input
            fi
        done
    fi



    # What was chosen: again, making difference between free choice and limited one
    echo "________________________________________________________________________________"
    echo ""
    echo "Analysis will be performed on :"
    if [[ "$4" == "fixed" ]]; then
        IFS=$'\n'
        declare -a nice_name_input=($( cat $8 ) )
        unset IFS
        # Print by sorted order, to keep the original order
        for each_input in $( echo $split_input | tr " " "\n" | sort -g ); do
            echo "    "${nice_name_input[$each_input]}
        done
    else
        # Print by sorted order
        for each_input in $input; do
            echo "    "$( cat $8 )" "$each_input
        done
    fi
}

# Confirming setup
confirm () {
setup_good=0
until [[ "$setup_good" -eq "1" ]]; do
    if [ -d $out_dir ]; then
        rm -r $out_dir
    fi
    setup $title "$intro" $table "$choice" $good_inputs $abbrev $out_dir $nice
    echo "Please press any key to confirm, or type [change] to change setup"
    read confirmation
    if [[ "$confirmation" != "change" ]]; then
        setup_good=1
        cd config
        rm -r temp/*
        cd ..
    fi
done
echo "________________________________________________________________________________"
echo "

"
}

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
table="config/temp/table.txt"
cat << "EOF" > $table
Type their corresponding number as single block of numbers (eg. 127), then press [ENTER]:

    1. Homo sapiens (long computational time!)
    2. Mus musculus (long computational time!)
    3. Caenorhabditis elegans
    4. Drosophila melanogaster
    5. Arabidopsis thaliana
    6. Saccharomyces cerevisiae
    7. Escherichia coli
EOF
choice=fixed
good_inputs=1234567
abbrev="config/temp/abbrev.txt"
cat << "EOF" > $abbrev
NA h_sapiens m_musculus c_elegans d_melanogaster a_thaliana s_cerevisiae e_coli
EOF
out_dir="config/species/"
nice="config/temp/nice.txt"
cat << "EOF" > $nice
NA
Homo sapiens
Mus musculus
Caenorhabditis elegans
Drosophila melanogaster
Arabidopsis thaliana
Saccharomyces cerevisiae
Escherichia coli
EOF

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



########### GENOMIC SIGNATURE TYPES ###########
title="config/temp/title.txt"
cat << "EOF" > $title
     _______                              __
    |     __|.-----.-----.-----.--------.|__|.----.
    |    |  ||  -__|     |  _  |        ||  ||  __|
    |_______||_____|__|__|_____|__|__|__||__||____|
            __                     __
    .-----.|__|.-----.-----.---.-.|  |_.--.--.----.-----.
    |__ --||  ||  _  |     |  _  ||   _|  |  |   _|  -__|
    |_____||__||___  |__|__|___._||____|_____|__| |_____|
               |_____|

EOF
intro="Please chose among the following list which type of genomic signature should be included in the analysis."
table="config/temp/table.txt"
cat << "EOF" > $table
Type their corresponding number as single block of numbers (eg. 12), then press [ENTER]:

    1. FCGR : k-mer frequencies using optimized computation (through the Chaos Game Representation
              (CGR) of a sequence)
    2. DFT : a lossless transformation of the CGR of a sequence, using Fourier transform
EOF
choice=fixed
good_inputs=12
abbrev="config/temp/abbrev.txt"
cat << "EOF" > $abbrev
NA FCGRs DFTs
EOF
out_dir="config/gs/"
nice="config/temp/nice.txt"
cat << "EOF" > $nice
NA
FCGR
DFT
EOF

confirm



########### FEATURES ###########
title="config/temp/title.txt"
cat << "EOF" > $title
 _______               __
|    ___|.-----.---.-.|  |_.--.--.----.-----.-----.
|    ___||  -__|  _  ||   _|  |  |   _|  -__|__ --|
|___|    |_____|___._||____|_____|__| |_____|_____|

EOF
intro="Please chose among the following list which features should be included in the analysis."
table="config/temp/table.txt"
cat << "EOF" > $table
Type their corresponding number as single block of numbers (eg. 12), then press [ENTER]:

    1. CDS : coding sequences
    2. LCR : low complexity regions (including transposable elements)
EOF
choice=fixed
good_inputs=12
abbrev="config/temp/abbrev.txt"
cat << "EOF" > $abbrev
NA CDS LCR
EOF
out_dir="config/features/"
nice="config/temp/nice.txt"
cat << "EOF" > $nice
NA
Coding sequences (CDS)
Low complexity regions (LCR)
EOF

confirm



########### K-MER SIZE ###########
# Checking if FCGR was chosen
if [ -f config/gs/FCGRs ]; then

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

fi



echo " ___ ___ _______ _______      ______ _______ ______ __  __ __
|   |   |       |   |   |    |   __ \       |      |  |/  |  |
 \     /|   -   |   |   |    |      <   -   |   ---|     <|__|
  |___| |_______|_______|    |___|__|_______|______|__|\__|__|
                                                              "

# Cleaning
cd config
rm -r temp
