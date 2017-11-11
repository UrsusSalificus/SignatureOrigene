#!/usr/bin/env bash

# Author : Titouan Laessle
# Copyright 2017 Titouan Laessle
# License : MIT

# Functions used in the setup part of the different analysis

# Reading the user input, creating the adequate config files and telling the user
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
    good_input=$( cat $5 | grep -o . )
else
    good_input=$( echo $5 | grep -o . )
fi

# Checking this is the right input
good_to_go=0
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
    if [[ -d $out_dir ]]; then
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