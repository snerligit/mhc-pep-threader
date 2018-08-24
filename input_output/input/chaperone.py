#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: August 22, 2018
#   Email: snerli@ucsc.edu
#

'''

CHAPERONE class contains all the necessary functionalities required to fetch
necessary chaperone sequences and their respective headers.

'''

# Load the Rosetta commands for use in the Python shell
from database.chaperone_sequences import chaperone_sequences
from input_output.input.header_readers import HEADER_READER

# inherit from HEADER_READER, a generic class that reads
# headers from a file
class CHAPERONE(HEADER_READER):

    # constructor
    def __init__(self, filename):
        super().__init__(filename)
        super().read_headers()

    # getter method
    def get_sequence(self, header):
        return chaperone_sequences[header]
