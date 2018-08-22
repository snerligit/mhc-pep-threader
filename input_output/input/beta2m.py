#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: August 22, 2018
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from database.beta2m_sequences import beta2m_sequences
from input_output.input.header_readers import HEADER_READER

class BETA2M(HEADER_READER):

    def __init__(self, filename):
        super().__init__(filename)
        super().read_headers()

    def get_sequence(self, header):
        return beta2m_sequences[header]
