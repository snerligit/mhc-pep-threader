#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 31, 2017
#   Email: snerli@ucsc.edu
#

# import other required libraries
import os
import sys
import math
import numpy
import argparse
import subprocess

# custom libraries
from movemap import MOVEMAP


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Perform threading of the template structure onto the target sequence')
    parser.add_argument("-pdb_file", help="provide template structure in PDB to perform threading", required='True')
    parser.add_argument("-groove_distance", help="provide the distance from peptide", default=3.5)
    parser.add_argument("-pep_start_index", help="provide start index of a peptide in the protein chain", type=int, required='True')
    parser.add_argument("-movemap_file_name", help="provide movemap filename", default='default.movemap')

    args = parser.parse_args()

    pdb_file_ = args.pdb_file
    groove_distance_ = args.groove_distance
    pep_start_index_ = args.pep_start_index
    movemap_file_name_ = args.movemap_file_name

    movemap = MOVEMAP(pdb_file_, pep_start_index_, groove_distance_, movemap_file_name_)
    movemap.apply()
