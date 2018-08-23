#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 26, 2017
#   Email: snerli@ucsc.edu
#

'''

ALIGN class contains all the necessary functionalities required to perform sequence
alignment between the template and the target sequences.

'''

# additional bio libraries
# Refernce: Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. 
# (2011 October 11) Molecular systems biology 7 :539
from Bio.Align.Applications import ClustalOmegaCommandline

#custom libraries
from alignment.substitution_matrix import BLOSUM

# import other required libraries
import os
import sys
import subprocess

class ALIGN:

    # class members
    template_seq = None # amino acid sequence for a given template
    target_seq = None # amino acid sequence for the required target
    matrix_type = 62 # matrix type to use for scoring alignment
    clustal_input = "" # input filename for clustal omega
    clustal_output = "" # output filename for clustal omega

    # constructor
    def __init__(self, template_seq, target_seq, matrix_type = 62):
        self.template_seq = template_seq
        self.matrix_type = matrix_type
        self.target_seq = target_seq
        self.clustal_input = "clustal_default_input.fasta"
        self.clustal_output = self.clustal_input+"_clustal_output.fasta"

    # method to create clustal input file
    def init_clustal_input(self):
        writefilehandle = open(self.clustal_input, "w")
        writefilehandle.write(">template\n")
        writefilehandle.write(self.template_seq+"\n")
        writefilehandle.write(">target\n")
        writefilehandle.write(self.target_seq+"\n")
        writefilehandle.close()

    # method that calls externally installed clustal program
    # read the input sequences from clustal input file and write to the
    # clustal output file
    def clustal(self):
        self.init_clustal_input()
        cline = ClustalOmegaCommandline(infile=self.clustal_input, outfile=self.clustal_output, distmat_full=True,
                                        verbose=True, seqtype="Protein", outfmt="vienna", iterations=10, percentid=True, force=True)
        cline()

    # getter method
    def get_clustal_output_filename(self):
        return self.clustal_output
