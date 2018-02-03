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

class FASTA:

    filename = ""
    header = ""
    seq = ""

    def __init__(self, filename, header="", seq=""):
        self.filename = filename
        self.header = header
        self.seq = seq

    def get_header(self):
        return self.header

    def get_seq(self):
        return self.seq

    def read(self):
        fastaFileHandler = open(self.filename, "r")
        for line in fastaFileHandler:
            line = line.rstrip()
            if ">" in line:
                self.header = "".join(list(line)[1:])
            else:
                self.seq += line

    def write(self):
        writefile = open(self.filename, "w")
        writefile.write(">"+self.header)
        writefile.write("\n")
        writefile.write(self.seq)
        writefile.close()
