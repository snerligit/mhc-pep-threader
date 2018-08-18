#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 26, 2017
#   Email: snerli@ucsc.edu
#

# additional bio libraries

#from skbio.core.alignment import StripedSmithWaterman
from skbio.alignment import StripedSmithWaterman
from Bio.Align.Applications import ClustalOmegaCommandline

#custom libraries
from alignment.substitution_matrix import BLOSUM

# import other required libraries
import os
import sys
import subprocess

class ALIGN:

    aln = None
    template_seq = None
    target_seq = None
    matrix_type = 62
    clustal_input = ""
    clustal_output = ""
    input_fasta = False

    def __init__(self, template_seq, target_seq, input_fasta = None, matrix_type = 62):
        self.template_seq = template_seq
        self.matrix_type = matrix_type
        if input_fasta == None and target_seq != None:
            self.target_seq = target_seq
            self.clustal_input = "clustal_default_input.fasta"
        else:
            self.clustal_input = input_fasta
            self.input_fasta = True
        self.clustal_output = self.clustal_input+"_clustal_output.fasta"

    def get_alignment(self):
        return self.aln

    def get_clustal_output_filename(self):
        return self.clustal_output

    def get_aligned_query(self):
        return self.aln.aligned_query_sequence

    def get_aligned_target(self):
        return self.aln.aligned_target_sequence

    def align(self):

        query = StripedSmithWaterman(self.template_seq, protein=True, substitution_matrix=BLOSUM(self.matrix_type).get_matrix())
        self.aln = query(self.target_seq)

    def check_if_template_exists(self):
        if self.input_fasta == True:
            readfilehandle = open(self.clustal_input, "r")
            for line in readfilehandle:
                line = line.rstrip()
                if ">template" in line:
                    return True
            readfilehandle.close()

    def init_clustal_input(self):
        if self.input_fasta == False:
            writefilehandle = open(self.clustal_input, "w")
            writefilehandle.write(">template\n")
            writefilehandle.write(self.template_seq+"\n")
            writefilehandle.write(">target\n")
            writefilehandle.write(self.target_seq+"\n")
            writefilehandle.close()
        else:
            if not self.check_if_template_exists():
                writefilehandle = open(self.clustal_input, "a")
                writefilehandle.write(">template\n")
                writefilehandle.write(self.template_seq+"\n")
                writefilehandle.close()


    def clustal(self):
        self.init_clustal_input()
        cline = ClustalOmegaCommandline(infile=self.clustal_input, outfile=self.clustal_output, distmat_full=True,
                                        verbose=True, seqtype="Protein", outfmt="vienna", iterations=10, percentid=True, force=True)
        cline()
        #self.read_aln()
