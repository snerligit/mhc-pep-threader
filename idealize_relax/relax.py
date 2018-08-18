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

    def relax_pdb(self, pose):
        relax = FastRelax()
        relax.set_scorefxn(get_fa_scorefxn())
        relax.apply(pose)
        return pose

    def relax_pdb_with_movemap(self, pose, movemap):
        relax = FastRelax()
        relax.set_scorefxn(get_fa_scorefxn())
        relax.set_movemap(movemap)
        relax.apply(pose)
        return pose
