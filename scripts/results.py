#!/usr/bin/env python3


import CGR_functions as fn
import subprocess

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

########################################################################################################################
                                        ############ 1 ###############
##############################
# - Window size: 150kb
# - Distance: DFTs + power spectrum of CGR -> euclidean distance
# - Feature: % of CDS
##############################

# Window size:
window_in_kb = '150kb'
# Type of genomic signature used:
gs = 'DFTs'
# Type of distance computed:
dist_type = 'euclidean'
# Type of feature:
feat_type = 'CDS'
# Features information:
feat_info = '% of'
# Paths to the distance matrices:
distances = fn.extract_path('/'.join(['../files/distances',dist_type, '_'.join([window_in_kb, gs]), '']), '*')
distances.sort()
# Paths to the features:
features = fn.extract_path('/'.join(['../files/features', feat_type, window_in_kb, '']), '*')
features.sort()

# Get all the species abbreviation for the study
species = fn.get_species()
species.sort()

for each_species in range(len(species)):
    outfile_cor = '/'.join(['../files/results', '_'.join([window_in_kb, gs, dist_type]),
                            '_'.join([species[each_species], 'correlation'])])
    outfile_MDS = '/'.join(['../files/results', '_'.join([window_in_kb, gs, dist_type]),
                            '_'.join([species[each_species], 'MDS'])])
    subprocess.call(['python', 'plot_PCA_and_cor.py', distances[each_species], features[each_species],
                     ' '.join([feat_info, feat_type]), '/'.join([window_in_kb, gs, dist_type]),
                     outfile_cor, outfile_MDS])
########################################################################################################################
