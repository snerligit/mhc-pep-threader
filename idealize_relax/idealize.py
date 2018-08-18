#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 26, 2017
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

# import rosetta files
from pyrosetta.rosetta.protocols.idealize import *
from pyrosetta.rosetta.core.scoring import *


class IDEALIZE:

    def idealize_pdb(self, pose):
        ideal = IdealizeMover()
        ideal.apply(pose)
        return pose
