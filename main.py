#!/usr/bin/env python3

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: February 2, 2018
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

#custom libraries
from model import MODEL
from input_output.input.argparser import ARGPARSE
from database.HLA_sequences_180 import hla_sequences_180
from database.HLA_sequences_complex import hla_sequences

# import other required libraries
import os
import sys
import subprocess

'''

Main is the entry point of RosettaMHC

'''

# method to run the RosettaMHC protocol
def run():
    args = ARGPARSE()

    if args.is_list_mhcs() == True:
        if not args.is_no_trim_mhc_flag_set():
            for key,value in hla_sequences.items():
                print(key)
        else:
            for key,value in hla_sequences_180.items():
                print(key)
    else:
        # Load Rosetta database files
        init(options='')
        modeller = MODEL(args)
        modeller.model_mhc_for_each_peptide_beta2m_tcr_chaperone()


if __name__ == "__main__":
    run()
