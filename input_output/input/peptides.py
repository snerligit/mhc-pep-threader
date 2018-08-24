#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: August 22, 2018
#   Email: snerli@ucsc.edu
#

'''

PEPTIDE class contains all the necessary functionalities required to fetch
necessary peptide sequences and their respective headers.

'''

from thread.fasta import FASTA

class PEPTIDE:

    # class members
    filename = "" # input file name containing peptide sequences in fasta format
    pep_object = None # peptide object of FASTA type

    # constructor
    def __init__(self, filename):
        self.filename = filename
        self.pep_object = FASTA(self.filename)
        self.pep_object.read()

    # getter methods
    def get_headers(self):
        return self.pep_object.get_headers()

    def get_sequences(self):
        return self.pep_object.get_sequences()

    def get_sequence(self, header):
        return self.pep_object.get_sequence(header)
