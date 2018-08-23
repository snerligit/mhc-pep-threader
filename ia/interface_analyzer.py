#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: February 1, 2018
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

class INTERFACE:

    pose = None
    ia = None

    def __init__(self, pose):
        self.ia = None
        self.pose = pose

    def analyze(self):
        self.ia = InterfaceAnalyzerMover()
        #self.ia.apply(self.pose)
        self.ia.apply_const(self.pose)

    def get_dG(self):
        return self.ia.get_interface_dG()
