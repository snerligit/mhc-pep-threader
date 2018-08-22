#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: August 22, 2018
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell

from thread.fasta import FASTA

class PEPTIDE:

    filename = ""
    pep_object = None

    def __init__(self, filename):
        self.filename = filename
        self.pep_object = FASTA(self.filename)
        self.pep_object.read()

    def get_headers(self):
        return self.pep_object.get_headers()

    def get_sequences(self):
        return self.pep_object.get_sequences()

    def get_sequence(self, header):
        return self.pep_object.get_sequence(header)
