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
from relax import RELAX
from template import TEMPLATE
from fasta import FASTA
from grishin import GRISHIN
from align import ALIGN
from movemap import MOVEMAP
from job_distributor import JOB_DISTRIBUTOR
from chain_split import CHAIN_SPLIT
from interface_analyzer import INTERFACE
from silent import SILENT
from thread import THREAD
from pre_thread import PRE_THREADING

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

'''
    def treatment_post_homology_modeling(self):

        threaded_pdb_name = self.template.get_template_path()+self.target_file_name+".pdb"
        self.threaded.dump_pdb(threaded_pdb_name)

        pep_start_index = find_pep_start_index()
        if self.movemap == None:
            movemap = MOVEMAP(threaded_pdb_name, self.pep_start_index, self.groove_distance, self.target_file_name+".movemap")
        else:
            movemap = MOVEMAP(threaded_pdb_name, self.pep_start_index, self.groove_distance, self.movemap, True)
        movemap_rosetta_obj = movemap.get_movemap()

        job_dist = JOB_DISTRIBUTOR(self.target_file_name, self.nstruct)

        relaxed_pose = self.threaded
        job_dist.get_dist_obj().native_pose = self.template.get_pose()

        while not job_dist.get_dist_obj().job_complete:
            relaxed_threaded_pose = RELAX(relaxed_pose)
            relaxed_threaded_pose.relax_pdb_with_movemap(movemap_rosetta_obj)
            job_dist.get_dist_obj().output_decoy(relaxed_pose)

            ia = INTERFACE(relaxed_pose)
            ia.analyze()
            print("Interface energy: ", ia.get_dG())

            relaxed_pose.assign(self.threaded)

'''
