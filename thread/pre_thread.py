#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: December 1, 2017
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

# import rosetta files
from pyrosetta.rosetta.protocols.comparative_modeling import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.core.sequence import *

#custom libraries
from idealize_relax.relax import RELAX
from thread.template import TEMPLATE
from thread.fasta import FASTA
from alignment.align import ALIGN
from alignment.grishin import GRISHIN

# import other required libraries
import os
import sys


'''

PRE_THREADING class contains all the necessary functionalities required to
perform operations before the threading process.

'''

class PRE_THREADING:

    # constructor
    template = None # store the template structure
    alignment = None # alignment object to store the alignment of template and target
    grishin_file_name = "" # grishin file name
    complex_header = "" # target complex header name
    complex_sequence = "" # target complex sequence
    pep_length = 0 # peptide length
    args = None # input arguments

    # constructor
    def __init__(self, template, complex_sequence, complex_header, args):
        self.template = template
        self.complex_sequence = complex_sequence
        self.complex_header = complex_header
        self.args = args
        self.pep_length = len(list(self.complex_header.split("_")[2]))
        self.align_template_target_sequences()

    # method to return target file name
    def get_target_file_name(self):
        return self.complex_header+"_on_"+self.template.get_stripped_name()

    # method to return target sequence
    def get_target_sequence(self):
        return self.complex_sequence

    # method to return peptide start index in the sequence
    def get_pep_start_index(self, sequence):
        index = 0
        try:
            index = sequence.rindex(self.peptide_sequence)
            return index
        except ValueError:
            print("Alignment failed. Cannot find peptide sequence")
            exit(1)

    # method to get the grishin file name
    def get_grishin_file_name(self):
        return self.grishin_file_name

    # method to return the template
    def get_template(self):
        return self.template

    # method to return the mhc header
    def get_mhc_header(self):
        return self.complex_header.split("_")[0]

    # method to return the peptide length
    def get_pep_length(self):
        return self.pep_length

    # method to check if grishin file exists
    def check_if_grishin_file_exists(self, filename):
        for f in self.grishin_file_name:
            if f == filename:
                return True

    # method to create a grishin file
    def create_grishin(self, aligned_target, aligned_query, is_new = False):
        # template_seq = query_sequence , target_seq =  target_sequence
        grishin = GRISHIN(self.get_target_file_name(), self.complex_header.split("_")[0], self.template.get_name(),
                            aligned_target, aligned_query)

        if not self.check_if_grishin_file_exists(grishin.get_file_name()):
            self.grishin_file_name = grishin.get_file_name()
            grishin.write()

    # method to align template and target sequences using clustal
    def align_template_target_sequences(self):
        key = self.complex_header.split("_")[0]
        self.alignment = ALIGN(self.template.get_sequence(), self.complex_header, self.complex_sequence, self.args.get_clustal_path())
        self.alignment.clustal()
        alignment = FASTA(self.alignment.get_clustal_output_filename())
        alignment.read()
        aln_seqs = alignment.get_sequences()

        for key,value in aln_seqs.items():
            if key != 'template':
                self.create_grishin(value, aln_seqs['template'])
