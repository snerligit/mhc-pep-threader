#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: April 10, 2017
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

# import rosetta files


# additional bio libraries

#custom libraries
from idealize_relax.relax import RELAX
from thread.template import TEMPLATE
from thread.fasta import FASTA
from idealize_relax.movemap import MOVEMAP
from jd.job_distributor import JOB_DISTRIBUTOR
from ia.chain_split import CHAIN_SPLIT
from ia.interface_analyzer import INTERFACE
from input_output.output.silent import SILENT
from thread.thread import THREAD
from thread.pre_thread import PRE_THREADING

from database.HLA_sequences_180 import hla_sequences_180
from database.HLA_sequences_complex import hla_sequences

# import other required libraries
import os
import sys


class MODEL_HLA:

    peptide_header = ""
    peptide_sequence = ""
    #beta2m = "MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
    tcr = ""
    mhc_fasta = None

    args = None

    def __init__(self, args):
        self.args = args
        if self.args.is_no_trim_mhc_flag_set():
            self.tcr = ""
        else:
            self.tcr = "KEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFWYRQYSGKSPELIMFIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSDSATYLCAVNFGGGKLIFGQGTELSVKPNIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSIAGITQAPTSQILAAGRRMTLRCTQDMRHNAMYWYRQDLGLGLRLIHYSNTAGTTGKGEVPDGYSVSRANTDDFPLTLASAVPSQTSVYFCASSLSFGTEAFFGQGTRLTVVEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQDPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD"

    def model_hlas_for_each_peptide(self):

        pep_file_handler = open(self.args.get_peptides(), "r")
        for line in pep_file_handler:
            line = line.rstrip()
            if ">" in line:
                self.peptide_header = line[1:]
            else:
                self.peptide_sequence = line
                self.apply()
        pep_file_handler.close()

    def apply(self):

        self.read_hla_sequences()
        pre_thread = PRE_THREADING(self.get_template(), self.mhc_fasta,
                            self.peptide_header, self.peptide_sequence,
                            self.args)
        threader = THREAD(pre_thread)
        threader.apply()
        #self.treatment_post_homology_modeling()

    def get_template(self):
        template_ = TEMPLATE(self.args.get_template_pdb())
        template_.treat_template_structure(self.args.get_mhc_chain(),
                                            self.args.get_peptide_chain(),
                                            self.args.is_no_trim_mhc_flag_set())
        return template_

    def read_from_map_and_write_to_fasta(self, fasta_file_handler, datamap):
        for key, value in datamap.items():
            fasta_file_handler.write(">"+key+"\n")
            fasta_file_handler.write(value+self.args.get_beta2m()+self.peptide_sequence+self.tcr+"\n")

    def read_key_from_map_and_write_to_fasta(self, fasta_file_handler, datamap, key):
        fasta_file_handler.write(">"+key+"\n")
        fasta_file_handler.write(datamap[key]+self.args.get_beta2m()+self.peptide_sequence+self.tcr+"\n")

    def read_hla_sequences(self):
        if self.args.get_target_fasta() == None:
            hla_file_handler = open(self.args.get_mhcs(), "r")
            self.mhc_fasta = self.args.get_mhcs()+".fasta"
            fasta_file_handler = open(self.mhc_fasta, "w")
            for line in hla_file_handler:
                line = line.rstrip()
                if line.lower() == "all":
                    if self.args.is_no_trim_mhc_flag_set():
                        self.read_from_map_and_write_to_fasta(fasta_file_handler, hla_sequences_180)
                    else:
                        self.read_from_map_and_write_to_fasta(fasta_file_handler, hla_sequences)
                else:
                    if self.args.is_no_trim_mhc_flag_set():
                        self.read_key_from_map_and_write_to_fasta(fasta_file_handler, hla_sequences_180, line)
                    else:
                        self.read_key_from_map_and_write_to_fasta(fasta_file_handler, hla_sequences, line)
            hla_file_handler.close()
            fasta_file_handler.close()
