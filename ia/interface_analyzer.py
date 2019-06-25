#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: February 1, 2018
#   Email: snerli@ucsc.edu
#

'''

INTERFACE class contains all the necessary functionalities required to calculate
interface energies between two chains in a pose.

'''

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

class INTERFACE:

    # class members
    pose = None # pose for which we want to calculte binding energies
    ia = None # InterfaceAnalyzerMover object that provides all the interface related interfaces

    # constructor
    def __init__(self, pose):
        self.ia = None
        self.pose = pose

    # method to apply the mover and populate necessary interface energies
    def analyze(self):
        self.ia = InterfaceAnalyzerMover()
        # we could have used self.ia.apply(self.pose), but the apply() function in turn calls
        # apply_const(), so we just call that directly
        self.ia.apply_const(self.pose)

    # getter method
    def get_dG(self):
        # return the interface score (or binding energy)
        return self.ia.get_interface_dG()
