#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 26, 2017
#   Email: snerli@ucsc.edu
#

# import other required libraries
import os
import sys
import math
import numpy
import argparse
import subprocess

# import biopython libraries
#
# Hamelryck, T., Manderick, B. (2003) PDB parser and structure class implemented in Python. Bioinformatics 19: 2308-2310
from Bio.PDB import *
#

class MyPDB:

    pdb = None
    filename = ""

    def __init__(self, filename):
        self.filename = filename

    def get_pdb(self):
        return self.pdb

    def read_chain_A(self):
        parser = PDBParser()
        structure = parser.get_structure(self.filename, self.filename)
        for model in structure:
            for chain in model:
                self.pdb = chain
                return
