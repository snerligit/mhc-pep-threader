#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 26, 2017
#   Email: snerli@ucsc.edu
#

# import biopython libraries
#
# Hamelryck, T., Manderick, B. (2003) PDB parser and structure class implemented in Python. Bioinformatics 19: 2308-2310
from Bio.PDB import *
#

'''

MyPDB class contains all the necessary functionalities required to
return the pdb( or required chain) given the filename.

'''

class MyPDB:

    # class members
    pdb = None # store the pdb read from PDB parser
    filename = "" # pdb filename

    # constructor
    def __init__(self, filename):
        self.filename = filename

    # getter method
    def get_pdb(self):
        return self.pdb

    # read only chain A of the input pdb
    def read_chain_A(self):
        parser = PDBParser()
        structure = parser.get_structure(self.filename, self.filename)
        for model in structure:
            for chain in model:
                self.pdb = chain
                return
