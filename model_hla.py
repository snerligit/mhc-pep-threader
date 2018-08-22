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
from thread.thread import THREAD
from thread.pre_thread import PRE_THREADING

from database.HLA_sequences_180 import hla_sequences_180
from database.HLA_sequences_complex import hla_sequences
from database.beta2m_sequences import beta2m_sequences
from database.tcr_sequences import tcr_sequences
from database.chaperone_sequences import chaperone_sequences

# import other required libraries
import os
import sys
import imp


class MODEL_HLA:

    peptide_headers = []
    peptide_sequences = []
    complex_headers = []
    complex_sequences = []
    template = None
    args = None
    mpi_install = True

    def __init__(self, args):
        self.args = args

    def model_hlas_for_each_peptide(self):

        pep_file_handler = open(self.args.get_peptides(), "r")
        for line in pep_file_handler:
            line = line.rstrip()
            if ">" in line:
                self.peptide_headers.append(line[1:])
            else:
                self.peptide_sequences.append(line)
                fasta_file_names = self.read_hla_sequences(line)
        pep_file_handler.close()

        njobs = len(self.complex_sequences)

        try:
            from jd.job_distributor import JOB_DISTRIBUTOR
            job_dist = JOB_DISTRIBUTOR()
            ''' #get this to work yet
            template_name = job_dist.perform_single_operation(self.get_template)
            if self.template == None and isinstance(template_name, str):
                self.template = TEMPLATE(template_name).get_pose_from_pdb()
                job_dist.apply(njobs, self.apply)
            elif self.template != None:
                job_dist.apply(njobs, self.apply)
            '''
            self.get_template()
            job_dist.apply(njobs, self.apply)
        except ImportError:
            self.mpi_install = False

        if not self.mpi_install:
            self.get_template()
            for job_id in range(njobs):
                self.apply(job_id)
                job_id += 1

    def apply(self, i):

        pre_thread = None
        pre_thread = PRE_THREADING(self.template, self.complex_sequences[i],
                            self.complex_headers[i], self.args)
        threader = THREAD(pre_thread)
        threader.apply()
        #self.treatment_post_homology_modeling()

    def get_template(self):
        self.template = TEMPLATE(self.args.get_template_pdb())
        self.template.treat_template_structure(self.args.get_mhc_chain(),
                                            self.args.get_peptide_chain(),
                                            self.args.is_no_trim_mhc_flag_set(),
                                            self.args.get_mhc_trim_length(),
                                            self.args.get_idealize_relax())
        return self.template.get_pdb()

    def write_to_fasta(self, key, value, peptide_sequence):
        fasta_file_name = key+"_"+peptide_sequence+".fasta"
        fasta_file_handler = open(fasta_file_name, "w")
        fasta_file_handler.write(">"+key+"\n")
        fasta_file_handler.write(value+"\n")
        fasta_file_handler.close()
        self.complex_headers.append(key+"_"+peptide_sequence)
        self.complex_sequences.append(value)
        return fasta_file_name

    def read_from_map_and_write_to_fasta(self, datamap, peptide_sequence, key = None):

        fasta_file_name = []

        if key == None:
            for key, value in datamap.items():
                text = value+beta2m_sequences[self.args.get_beta2m()]+peptide_sequence+tcr_sequences[self.args.get_tcr()]+chaperone_sequences[self.args.get_chaperone()]
                fasta_file_name.append(self.write_to_fasta(key, text, peptide_sequence))
        else:
            text = datamap[key]+beta2m_sequences[self.args.get_beta2m()]+peptide_sequence+tcr_sequences[self.args.get_tcr()]+chaperone_sequences[self.args.get_chaperone()]
            fasta_file_name.append(self.write_to_fasta(key, text, peptide_sequence))

        return fasta_file_name

    def read_hla_sequences(self, peptide_sequence):
        fasta_file_names = []
        if self.args.get_target_fasta() == None:
            hla_file_handler = open(self.args.get_mhcs(), "r")
            for line in hla_file_handler:
                line = line.rstrip()
                key = line
                if line.lower() == "all":
                    key = None
                if self.args.is_no_trim_mhc_flag_set():
                    fasta_file_names += (self.read_from_map_and_write_to_fasta(hla_sequences_180, peptide_sequence, key))
                else:
                    fasta_file_names += (self.read_from_map_and_write_to_fasta(hla_sequences, peptide_sequence, key))
            hla_file_handler.close()
        return fasta_file_names
