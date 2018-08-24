#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: February 1, 2018
#   Email: snerli@ucsc.edu
#

'''

CHAIN_SPLIT class contains all the necessary functionalities required to divide
a pose into two separate chains based on the given residue number also referred to as cutpoint.

'''

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_moves import ChainSplitMover

class CHAIN_SPLIT:

    # class members
    pose = None # pose which needs to be split into two chains based on the residue number
    cutpoint = 1 # residue number that marks the split between two chains

    # constructor
    def __init__(self, pose, cutpoint):
        self.pose = pose
        self.cutpoint = cutpoint

    # method to perform chain split
    def cut(self):
        split = ChainSplitMover(self.cutpoint)
        split.apply(self.pose)

    # getter method
    def get_pose(self):
        return self.pose
