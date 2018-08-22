#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 26, 2017
#   Email: snerli@ucsc.edu
#

# additional bio libraries

# import other required libraries
import os
import sys
import subprocess
from collections import defaultdict

class FASTA:

    filename = ""
    sequences = defaultdict(dict)

    def __init__(self, filename):
        self.filename = filename
        self.sequences = defaultdict(dict)

    def get_headers(self):
        if not self.sequences:
            self.sequences['none'] = ''
        return self.sequences.keys()

    def get_sequences(self):
        return self.sequences

    def get_sequence(self, header):
        return self.sequences[header]

    def read(self):
        if self.filename != None:
            fasta_file_handler = open(self.filename, "r")
            for line in fasta_file_handler:
                line = line.rstrip()
                if ">" in line:
                    header = "".join(list(line)[1:])
                else:
                    self.sequences[header] = line
            fasta_file_handler.close()

    def write(self, header, text):
        fasta_file_handler = open(self.filename, "w")
        fasta_file_handler.write(">"+header+"\n")
        fasta_file_handler.write(text+"\n")
        fasta_file_handler.close()
