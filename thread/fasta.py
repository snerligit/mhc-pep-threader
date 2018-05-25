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

    def get_headers(self):
        return self.sequences.keys()

    def get_sequences(self):
        return self.sequences

    def get_sequence(self, header):
        return self.sequences[header]

    def read(self):
        fastaFileHandler = open(self.filename, "r")
        for line in fastaFileHandler:
            line = line.rstrip()
            if ">" in line:
                header = "".join(list(line)[1:])
            else:
                self.sequences[header] = line
