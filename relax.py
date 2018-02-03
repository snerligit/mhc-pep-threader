#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 26, 2017
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

# import rosetta files
from pyrosetta.rosetta.protocols.relax import *
from pyrosetta.rosetta.core.scoring import *


class RELAX:

    pose = None

    def __init__(self, pose):
        self.pose = pose

    def relax_pdb(self):
        relax = FastRelax()
        relax.set_scorefxn(get_fa_scorefxn())
        relax.apply(self.pose)
        return pose

    def relax_pdb_with_movemap(self, movemap):
        relax = FastRelax()
        relax.set_scorefxn(get_fa_scorefxn())
        relax.set_movemap(movemap)
        relax.apply(self.pose)
        return self.pose
