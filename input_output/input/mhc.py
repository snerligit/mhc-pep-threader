#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: August 22, 2018
#   Email: snerli@ucsc.edu
#

'''

MHC class contains all the necessary functionalities required to fetch
necessary mhc sequences and their respective headers.

'''

# Load the Rosetta commands for use in the Python shell
from database.HLA_sequences_180 import hla_sequences_180
from database.HLA_sequences_complex import hla_sequences
from input_output.input.header_readers import HEADER_READER

# inherit from HEADER_READER, a generic class that reads
# headers from a file
class MHC(HEADER_READER):

    # class member
    trim_mhc = True # to check if mhc sequence needs to be trimmed

    # constructor
    def __init__(self, filename, trim_mhc):
        super().__init__(filename)
        super().read_headers()
        self.trim_mhc = trim_mhc

    # getter methods
    # based on the trim flag, a full sequencce or a trimmed
    # sequence is returned
    def get_sequence(self, header):
        if self.trim_mhc == True:
            return hla_sequences_180[header]
        else:
            return hla_sequences[header]

    # a user can just specify all in the input mhc list file
    # if it is specified, a full list of headers is returned
    # otherwise, only a requested subset of headers is returned
    def get_headers(self):
        if super().get_headers()[0].lower != "all":
            return super().get_headers()
        else:
            if self.trim_mhc == True:
                return hla_sequences_180.keys()
            else:
                return hla_sequences.keys()
