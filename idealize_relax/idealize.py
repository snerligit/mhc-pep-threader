#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 26, 2017
#   Email: snerli@ucsc.edu
#

'''

IDEALIZE class contains all the necessary functionalities required to idealize
an input pose.

'''

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

# import rosetta files
from pyrosetta.rosetta.protocols.idealize import *

class IDEALIZE:

    # method to apply IdealizeMover for a given pose
    def idealize_pdb(self, pose):
        ideal = IdealizeMover()
        ideal.apply(pose)
        return pose
