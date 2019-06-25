#!/usr/bin/python

#       Sgourakis Lab
#   Author: Hailey Wallace
#   Date: June 19, 2019
#   Email: hmwalla2@nscu.edu

'''

BACKRUB class contains all the necessary functionalities required to backrub
an input pose.

'''

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

# import rosetta files
from pyrosetta.rosetta.protocols.backrub import *

class BACKRUB:

    # method to apply BackrubProtocol for a given pose
    def backrub_pdb(self, pose):
        br = BackrubProtocol()
        br.set_scorefunction(get_fa_scorefxn())
        br.apply(pose)
        return pose

    #method to perturb the backbone using Backrub
    #also perturbs those residues specified in the movemap

    def backrub_pdb_with_movemap(self, pose, movemap):
        br = BackrubProtocol()
        br.set_scorefunction(get_fa_scorefxn())
        br.set_movemap(movemap)
        br.apply(pose)
        return pose
