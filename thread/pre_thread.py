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

# additional bio libraries

#custom libraries
from idealize_relax.relax import RELAX
from thread.template import TEMPLATE
from thread.fasta import FASTA
from alignment.align import ALIGN
from alignment.grishin import GRISHIN

# import other required libraries
import os
import sys

class PRE_THREADING:

    template = None
    alignment = None
    grishin_file_name = ""
    complex_header = ""
    complex_sequence = ""
    #pep_start_index = 181
    #pep_start_index = 376
    pep_length = 0
    args = None

    def __init__(self, template, complex_sequence, complex_header, args):

        self.template = template
        self.complex_sequence = complex_sequence
        self.complex_header = complex_header
        self.args = args
        self.pep_length = len(list(self.complex_header.split("_")[1]))
        self.align_template_target_sequences()

    def get_target_file_name(self):
            fields = self.complex_header.split("_")
            mhc_header = fields[0]
            peptide_header = fields[1]
            return mhc_header+"_on_"+self.template.get_stripped_name()+"_with_"+peptide_header

    def get_target_sequence(self):
        
        return self.complex_sequence

    def get_pep_start_index(self, sequence):

        index = 0
        try:
            index = sequence.rindex(self.peptide_sequence)
            return index
        except ValueError:
            print("Retry with different alignment scheme")
            exit(1)

    def check_if_grishin_file_exists(self, filename):
        for f in self.grishin_file_name:
            if f == filename:
                return True

    def create_grishin(self, aligned_target, aligned_query, is_new = False):

        # template_seq = query_sequence , target_seq =  target_sequence
        grishin = GRISHIN(self.get_target_file_name(), self.complex_header.split("_")[0], self.template.get_name(),
                            aligned_target, aligned_query)

        if not self.check_if_grishin_file_exists(grishin.get_file_name()):
            self.grishin_file_name = grishin.get_file_name()
            grishin.write()

    def get_grishin_file_name(self):

        return self.grishin_file_name

    def get_template(self):

        return self.template

    def get_mhc_header(self):

        return self.complex_header.split("_")[0]

    def get_pep_length(self):

        return self.pep_length

    def align_template_target_sequences(self, clustal=True):

        key = self.complex_header.split("_")[0]
        if clustal == False:
            self.alignment = ALIGN(self.template.get_sequence(), self.complex_sequence)
            self.alignment.align()
            self.create_grishin(self.alignment.get_aligned_target(), self.alignment.get_aligned_query())

        else:
            self.alignment = ALIGN(self.template.get_sequence(), self.complex_sequence)
            self.alignment.clustal()
            alignment = FASTA(self.alignment.get_clustal_output_filename())
            alignment.read()
            aln_seqs = alignment.get_sequences()

            for key,value in aln_seqs.items():
                if key != 'template':
                    self.create_grishin(value, aln_seqs['template'])
