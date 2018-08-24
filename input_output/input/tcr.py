#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: August 22, 2018
#   Email: snerli@ucsc.edu
#

'''

TCR class contains all the necessary functionalities required to fetch
necessary tcr sequences and their respective headers.

'''

from thread.fasta import FASTA

class TCR:

    # class members
    filename = "" # input file name containing tcr sequences in fasta format
    tcr_object = None # tcr object of FASTA type

    # constructor
    def __init__(self, filename):
        self.filename = filename
        self.tcr_object = FASTA(self.filename)
        self.tcr_object.read()

    # getter methods
    def get_headers(self):
        return self.tcr_object.get_headers()

    def get_sequences(self):
        return self.tcr_object.get_sequences()

    def get_sequence(self, header):
        return self.tcr_object.get_sequence(header)
