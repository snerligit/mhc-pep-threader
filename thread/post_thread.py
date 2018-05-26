#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: May 25, 2017
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

class POST_THREADING:

    args = None

    def __init__(self, args):

        self.args = args

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
