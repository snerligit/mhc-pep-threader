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
from thread.post_thread import POST_THREADING

# import other required libraries
import os
import sys

'''

THREAD class contains all the necessary functionalities required to
perform target onto the template structure.

'''

class THREAD:

    # class members
    pre_threader = None # pre threader object

    # constructor
    def __init__(self, pre_thread):
        self.pre_threader = pre_thread

    # method that calls threader
    def apply(self):
        if self.pre_threader != None:
            self.do_partial_threading()

    # method to perform partial threading
    def do_partial_threading(self):
        # for each alignment combination in the grishin file.
        alignment_vec = read_aln("grishin", self.pre_threader.get_grishin_file_name())
        for vec in alignment_vec:
            new_pose = pose_from_sequence(self.pre_threader.get_target_sequence(), "fa_standard")
            thread = PartialThreadingMover(vec, self.pre_threader.get_template().get_pose())

            thread.apply(new_pose)

            threaded_pose = new_pose
            tag = self.pre_threader.get_target_file_name()

            post = POST_THREADING(self.pre_threader, threaded_pose, tag)
            post.treatment_post_homology_modeling()
