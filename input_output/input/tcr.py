#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: August 22, 2018
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell

from thread.fasta import FASTA

class TCR:

    filename = ""
    tcr_object = None

    def __init__(self, filename):
        self.filename = filename
        self.tcr_object = FASTA(self.filename)
        self.tcr_object.read()

    def get_headers(self):
        return self.tcr_object.get_headers()

    def get_sequences(self):
        return self.tcr_object.get_sequences()

    def get_sequence(self, header):
        return self.tcr_object.get_sequence(header)
