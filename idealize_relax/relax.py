#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 26, 2017
#   Email: snerli@ucsc.edu
#

'''

RELAX class contains all the necessary functionalities required to call
the appropriate minimizer for a given pose.

'''

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

# import rosetta files
from pyrosetta.rosetta.protocols.relax import *
from pyrosetta.rosetta.core.scoring import *

class RELAX:

    # method to perform simple minimization using FastRelax
    def relax_pdb(self, pose):
        relax = FastRelax()
        relax.set_scorefxn(get_fa_scorefxn())
        relax.apply(pose)
        return pose

    # method to perform simple minimization using FastRelax
    # also minimizes only those residues specified in the movemap
    def relax_pdb_with_movemap(self, pose, movemap):
        relax = FastRelax()
        relax.set_scorefxn(get_fa_scorefxn())
        relax.set_movemap(movemap)
        relax.apply(pose)
        return pose
