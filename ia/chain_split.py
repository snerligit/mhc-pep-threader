#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: February 1, 2018
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_moves import ChainSplitMover

class CHAIN_SPLIT:

    pose = None
    cutpoint = 1

    def __init__(self, pose, cutpoint):
        self.pose = pose
        self.cutpoint = cutpoint

    def cut(self):
        split = ChainSplitMover(self.cutpoint)
        split.apply(self.pose)

    def get_pose(self):
        return self.pose
