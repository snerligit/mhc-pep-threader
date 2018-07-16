#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: February 2, 2018
#   Email: snerli@ucsc.edu
#


# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

#custom libraries
from model_hla import MODEL_HLA
from input_output.input.argparser import ARGPARSE
#from database.HLA_sequences_180 import hla_sequences_180
from database.HLA_sequences import hla_sequences
# added b2m sequences - ACM
from database.beta2m_sequences import beta2m_sequences

# import other required libraries
import os
import sys
import subprocess

def run():
    args = ARGPARSE()

    if args.is_list_mhcs() == True:
        if args.is_no_trim_mhc_flag_set():
            for key,value in hla_sequences.items():
                print(key)
        else:
            # changed to hla_sequences_180 - ACM
            for key,value in hla_sequences.items():
                print(key)
    else:
        # Load Rosetta database files
        init()
        modeller = MODEL_HLA(args)
        modeller.model_hlas_for_each_peptide()


if __name__ == "__main__":
    run()
