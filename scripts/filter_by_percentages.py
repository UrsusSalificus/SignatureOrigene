"""This script will extract the overall percentages of each factor
"""
import sys
import os.path

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

# Factor percentages for this species:
factor_percentages = str(sys.argv[1])
# Species name
species = '_'.join(str(factor_percentages.split('/')[-1]).split('_')[:2])
# Follow up:
follow_up = str(sys.argv[2])


###
# Check if parent directory is present, if not create it
###
def checking_parent(file_path):
    # We don't need the file name, so will take everything but the last part
    parent_directories = '/'.join(file_path.split('/')[0:(len(file_path.split('/')) - 1)])
    # As we uses parallel, we ended up we one thread doesn't seeing the directory, attempting
    # creating it, while another just did the same -> error "The file already exist", and stopped everything...
    try:
        if not os.path.exists(parent_directories):
            os.makedirs(parent_directories)
    except:
        pass


# This list will only keep the factors present above a certain threshold
factor_kept = list()

with open(factor_percentages, 'r') as percs:
    for each_line in percs:
        line = each_line.strip().split('\t')
        if float(line[2]) > 0.1 and int(line[4]) >= 3:
            factor_kept.append(line[1])

# Write these new factors in the config file
new_config_directory = 'config/new_factors/' + species + '/'

# Writing the factors kept
checking_parent(new_config_directory)
for each_factor in factor_kept:
    with open(new_config_directory + each_factor, 'w') as outfile:
        outfile.write('')

# Follow the progression of the analysis
checking_parent(follow_up)
with open(follow_up, 'w') as file:
    file.write('')