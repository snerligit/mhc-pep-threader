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
from idealize_relax.movemap import MOVEMAP

# import other required libraries
import os
import sys


class THREAD:

    pre_threader = None

    def __init__(self, pre_thread):

        self.pre_threader = pre_thread

    def apply(self):

        self.do_partial_threading()

    def do_partial_threading(self):

        alignment_vec = read_aln("grishin", self.pre_threader.get_grishin_file_name())
        ctr = 0
        for vec in alignment_vec:
            new_pose = pose_from_sequence(self.pre_threader.get_target().get_sequence(self.pre_threader.get_mhc_header(ctr)), "fa_standard")
            thread = PartialThreadingMover(vec, self.pre_threader.get_template().get_pose())
            thread.apply(new_pose)

            threaded = new_pose
            tag = self.pre_threader.get_target_file_name(self.pre_threader.get_mhc_header(ctr))
            threaded.dump_pdb(tag+".pdb")
            movemap = MOVEMAP(tag+".pdb", self.pre_threader.args.get_pep_start_index(),
                                self.pre_threader.get_pep_length(), self.pre_threader.args.get_groove_distance(),
                                tag+".movemap", True)
            movemap.apply()
            ctr += 1
