#!/usr/bin/env python3

"""Plot both the PCA and xy plot of distance matrices and features
"""
import CGR_functions as fn
import subprocess
import sys

__author__ = "Titouan Laessle"
__copyright__ = "Copyright 2017 Titouan Laessle"
__license__ = "MIT"

fn.checking_parent(sys.argv[5])

subprocess.call(['Rscript', 'cor_dist_feature.R', sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]])
subprocess.call(['Rscript', 'MDS_dist_feature.R', sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[6]])