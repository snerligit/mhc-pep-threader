#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: February 2, 2018
#   Email: snerli@ucsc.edu
#


# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

#custom libraries
from thread import THREAD

# import other required libraries
import os
import sys
import subprocess
import argparse

# Load Rosetta datbase files
init()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Perform threading of the template structure onto the target sequence')
    parser.add_argument("-template_pdb", help="provide template structure in PDB to perform threading", required='True')
    parser.add_argument("-fasta", help="provide FASTA format file with single sequence of the target being modeled", required='True')
    parser.add_argument("-movemap", help="provide a movemap file for relaxing the target")
    parser.add_argument("-groove_distance", help="provide distance to select nearest groove residues from the peptide", type=float, default=3.5)
    parser.add_argument("-pep_start_index", help="provide start index of a peptide in the protein chain", type=int, required='True')
    parser.add_argument("-idealize_relax", help="idealize and relax template structure before threading", default=False, type=bool)
    parser.add_argument("-nstruct", help="number of times a threaded structure should be relaxed", type=int, default=1)
    parser.add_argument("-target_file_name", help="provide a file name for the threaded structure")

    args = parser.parse_args()

    template_pdb_ = args.template_pdb
    fasta_ = args.fasta
    idealize_relax_ = args.idealize_relax
    nstruct_ = args.nstruct
    movemap_ = args.movemap
    groove_distance_ = args.groove_distance
    pep_start_index_ = args.pep_start_index
    target_file_name = args.target_file_name

    threader = THREAD(template_pdb_, fasta_, pep_start_index_, groove_distance_, nstruct_, movemap_, target_file_name)
    threader.apply()
