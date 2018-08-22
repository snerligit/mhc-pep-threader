#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: August 22, 2018
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from database.HLA_sequences_180 import hla_sequences_180
from database.HLA_sequences_complex import hla_sequences
from input_output.input.header_readers import HEADER_READER

class MHC(HEADER_READER):

    trim_mhc = True

    def __init__(self, filename, trim_mhc):
        super().__init__(filename)
        super().read_headers()
        self.trim_mhc = trim_mhc

    def get_sequence(self, header):
        if self.trim_mhc == True:
            return hla_sequences_180[header]
        else:
            return hla_sequences[header]

    def get_headers(self):
        if super().get_headers()[0].lower != "all":
            return super().get_headers()
        else:
            if self.trim_mhc == True:
                return hla_sequences_180.keys()
            else:
                return hla_sequences.keys()
