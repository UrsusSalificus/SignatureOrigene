#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# Will setup the whole analysis through user commands
# ASCII arty letters created through    http://patorjk.com/software/taag/#p=display&f=Chunky&t=Example

rm -r config/*
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
echo "________________________________________________________________________________________________________________"
echo ""
echo Welcome to the setup of the analysis
echo "________________________________________________________________________________________________________________"
echo ""
#read -n 1 -s -r -p "Please, press any key to continue"
echo ""
echo "________________________________________________________________________________________________________________"
echo "


"

########### DEFINITION OF THE CHOOSING FUNCTION ###########
setup () {
    cat $1
    echo $2
    echo "                  ### To do so ###"
    echo "Type their corresponding number (eg. 12), then press [ENTER]"
    cat $3
    read input

    # Checking not empty
    while [[ "$input" == "" ]]; do
      echo "Please, enter at least one input."
      read input
    done

    split_input=$( echo $input | grep -o . )

    # Checking this is the right input
    good_to_go=0
    good_input=$( echo $4 | grep -o . )
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
    declare -a all_input=($( cat $5 ))
    for each_input in $split_input; do
        if [ ! -f $6${all_input[$each_input]} ]; then
            mkdir -p $6 && touch $6${all_input[$each_input]}
        fi
    done

    # What was chosen:
    IFS=$'\n'
    declare -a nice_name_input=($( cat $7 ) )
    unset IFS
    echo "________________________________________________________________________________________________________________"
    echo ""
    echo "Analysis will be performed on :"
    # Print by sorted order, to keep the original order
    for each_input in $( echo $split_input | tr " " "\n" | sort -g ); do
        echo ${nice_name_input[$each_input]}
    done
}

# Confirming setup
confirm () {
setup_good=0
until [[ "$setup_good" -eq "1" ]]; do
    if [ -d $out_dir ]; then
        rm -r $out_dir
    fi
    setup $title "$intro" $table $good_inputs $abbrev $out_dir $nice
    echo "Please press any key to confirm, or type [change] to change setup"
    read confirmation
    if [[ "$confirmation" != "change" ]]; then
        setup_good=1
        cd config
        rm -r temp/*
        cd ..
    fi
done
echo "________________________________________________________________________________________________________________"
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
    1. Homo sapiens (long computational time!)
    2. Mus musculus (long computational time!)
    3. Caenorhabditis elegans
    4. Drosophila melanogaster
    5. Arabidopsis thaliana
    6. Saccharomyces cerevisiae
    7. Escherichia coli
EOF
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
    1. FCGR : k-mer frequencies using optimized computation (through the Chaos Game Representation
              (CGR) of a sequence)
    2. DFT : a lossless transformation of the CGR of a sequence, using Fourier transform
EOF
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


# Cleaning
cd config
rm -r temp
