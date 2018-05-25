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
from relax import RELAX
from template import TEMPLATE
from fasta import FASTA
from grishin import GRISHIN
from align import ALIGN

# import other required libraries
import os
import sys

class PRE_THREADING:

    target = None
    template = None
    alignment = None
    grishin_file_name = None
    mhc_fasta_file = ""
    peptide_sequence = ""
    peptide_header = ""
    mhc_headers = []
    #pep_start_index = 181
    #pep_start_index = 376
    pep_length = 0
    args = None

    def __init__(self, template, mhc_fasta_file, peptide_header, peptide_sequence, args):

        self.template = template
        self.mhc_fasta_file = mhc_fasta_file
        self.peptide_sequence = peptide_sequence
        self.peptide_header = peptide_header
        self.args = args
        self.pep_length = len(list(peptide_sequence))
        self.align_template_target_sequences()

    def get_target_file_name(self, header="HLAs"):

            return header+"_on_"+self.template.get_stripped_name()+"_with_"+self.peptide_header

    def get_pep_start_index(self, sequence):

        index = 0
        try:
            index = sequence.rindex(self.peptide_sequence)
            return index
        except ValueError:
            print("Retry with different alignment scheme")
            exit(1)

    def create_grishin(self, key, aligned_target, aligned_query):

        # template_seq = query_sequence , target_seq =  target_sequence
        grishin = GRISHIN(self.get_target_file_name(), key, self.template.get_name(),
                            aligned_target, aligned_query)
        self.grishin_file_name = grishin.get_file_name()
        grishin.write()

    def get_grishin_file_name(self):

        return self.grishin_file_name

    def get_template(self):

        return self.template

    def get_target(self):

        return self.target

    def get_mhc_header(self, key):

        return self.mhc_headers[key]

    def get_pep_length(self):

        return self.pep_length

    def align_template_target_sequences(self, clustal=True):

        self.target = FASTA(self.mhc_fasta_file)
        if clustal == False:
            self.target.read()
            target_seqs = self.target.get_sequences()

            for key,value in target_seqs.items():
                self.alignment = ALIGN(self.template.get_sequence(), self.target.get_sequence(key))
                self.alignment.align()
                self.create_grishin(key, self.alignment.get_aligned_target(), self.alignment.get_aligned_query())
                self.mhc_headers.append(key)

        else:
            print(self.args.get_target_fasta())
            self.alignment = ALIGN(self.template.get_sequence(), None, self.mhc_fasta_file)
            self.alignment.clustal()
            alignment = FASTA(self.alignment.get_clustal_output_filename())
            alignment.read()
            aln_seqs = alignment.get_sequences()

            for key,value in aln_seqs.items():
                if key != 'template':
                    self.create_grishin(key, value, aln_seqs['template'])
                    self.mhc_headers.append(key)
