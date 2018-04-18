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

from database.HLA_sequences_180 import hla_sequences

# import other required libraries
import os
import sys


class MODEL_HLA:

    target = None
    template = None
    alignment = None
    grishin_file_name = None
    threaded = None
    template_pdb = None
    fasta = ""
    groove_distance = 3.5
    nstruct = 1
    movemap = None
    peptide = ""
    pep_header = ""
    #beta2m = "MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
    beta2m = ""

    def __init__(self, template_pdb, template, fasta, groove_distance, nstruct, movemap, pepite_header, peptide, hlas, beta2m):

        self.template_pdb = template_pdb
        self.template = template
        self.fasta = fasta
        self.groove_distance = groove_distance
        self.nstruct = nstruct
        self.movemap = movemap
        self.pep_header = pepite_header
        self.peptide = peptide
        self.hlas = hlas
        self.beta2m = beta2m

    def apply(self):

        self.read_hla_sequences()
        threader = THREAD(self.template, self.fasta, self.nstruct, self.pep_header, self.peptide)
        threader.apply()
        #self.treatment_post_homology_modeling()

    def read_hla_sequences(self):
        if self.fasta == None:
            hla_file_handler = open(self.hlas, "r")
            self.fasta = self.hlas+".fasta"
            fasta_file_handler = open(self.fasta, "w")
            for line in hla_file_handler:
                line = line.rstrip()
                if line.lower() == "all":
                    for key, value in hla_sequences:
                        fasta_file_handler.write(">"+key+"\n")
                        fasta_file_handler.write(value+self.beta2m+self.peptide+"\n")
                else:
                    fasta_file_handler.write(">"+line+"\n")
                    fasta_file_handler.write(hla_sequences[line]+self.beta2m+self.peptide+"\n")
                    print(hla_sequences[line]+self.peptide+"\n")
            hla_file_handler.close()
            fasta_file_handler.close()

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
