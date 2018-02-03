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
from pyrosetta.rosetta.protocols.jobdist import *

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

# import other required libraries
import os
import sys


class THREAD:

    target = None
    target_file_name = "default"
    template = None
    alignment = None
    grishin_file_name = None
    threaded = None
    template_pdb = None
    fasta = ""
    pep_start_index = 1
    groove_distance = 3.5
    nstruct = 1
    movemap = None

    def __init__(self, template_pdb, fasta, pep_start_index, groove_distance, nstruct, movemap, target_file_name):

        self.template_pdb = template_pdb
        self.fasta = fasta
        self.pep_start_index = pep_start_index
        self.groove_distance = groove_distance
        self.nstruct = nstruct
        self.movemap = movemap
        self.target_file_name = target_file_name

    def apply(self):

        self.template = TEMPLATE(self.template_pdb)
        self.template.treat_template_structure()
        self.align_template_target_sequences()
        self.do_homology_modeling()
        self.treatment_post_homology_modeling()

    def align_template_target_sequences(self):

        self.target = FASTA(self.fasta)
        self.target.read()
        self.set_target_file_name()

        self.alignment = ALIGN(self.template.get_seq(), self.target.get_seq())
        self.alignment.align()

        # template_seq = query_sequence , target_seq =  target_sequence
        grishin = GRISHIN(self.template.get_template_path()+self.target_file_name, self.template.get_name(),
                self.target.get_header(), self.alignment.get_aligned_query(),
                self.alignment.get_aligned_target())
        grishin.write()

        self.grishin_file_name = grishin.get_file_name()

    def set_target_file_name(self):

        if self.target_file_name == None:
            self.target_file_name = self.target.get_header()+"_on_"+self.template.get_stripped_name()

    def do_homology_modeling(self):

        alignment_vec = read_aln("grishin", self.grishin_file_name)
        new_pose = pose_from_sequence(self.target.get_seq(), "fa_standard")
        for vec in alignment_vec:
            print(vec)
            thread = PartialThreadingMover(vec, self.template.get_pose())
            thread.apply(new_pose)

        chain = CHAIN_SPLIT(new_pose, self.pep_start_index)
        chain.cut()
        self.threaded = chain.get_pose()

    def treatment_post_homology_modeling(self):

        threaded_pdb_name = self.template.get_template_path()+self.target_file_name+".pdb"
        self.threaded.dump_pdb(threaded_pdb_name)

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
            print("Centoid dG: ", ia.get_dG())

            relaxed_pose.assign(self.threaded)
