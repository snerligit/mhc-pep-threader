#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 26, 2017
#   Email: snerli@ucsc.edu
#

# additional bio libraries

#from skbio.core.alignment import StripedSmithWaterman
from skbio.alignment import StripedSmithWaterman

#custom libraries
from substitution_matrix import BLOSUM

# import other required libraries
import os
import sys
import subprocess

class ALIGN:

    aln = None
    template_seq = None
    target_seq = None
    matrix_type = 62

    def __init__(self, template_seq, target_seq, matrix_type = 62):
        self.template_seq = template_seq
        self.target_seq = target_seq
        self.matrix_type = matrix_type

    def get_alignment(self):
        return self.aln

    def get_aligned_query(self):
        return self.aln.aligned_query_sequence

    def get_aligned_target(self):
        return self.aln.aligned_target_sequence

    def align(self):

        query = StripedSmithWaterman(self.template_seq, protein=True, substitution_matrix=BLOSUM(self.matrix_type).get_matrix())
        self.aln = query(self.target_seq)
