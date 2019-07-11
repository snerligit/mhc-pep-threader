#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: July 9, 2018
#   Email: snerli@ucsc.edu
#

'''

RMSD class contains all the necessary functionalities required to compute RMSD values
of the threaded/relaxed models to the input native structure.

'''

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *
from pyrosetta.rosetta.core.scoring import *

class RMSD:

    # class members
    pose = None # pose whose rmsd values are calculated
    native_pose = None # native pose
    residues_of_interest = [] # residues used to compute rmsds over

    # constructor
    def __init__(self, pose, native_pose, residues_of_interest = None):
        self.pose = pose
        self.native_pose = native_pose
        self.residues_of_interest = residues_of_interest

    def CA_rmsd(self):
        rms_ca = pyrosetta.rosetta.core.scoring.CA_rmsd(self.native_pose, self.pose)
        return rms_ca

    def CA_rmsd_between_indices(self, start, end):
        rms_ca = pyrosetta.rosetta.core.scoring.CA_rmsd(self.native_pose, self.pose, start, end)
        return rms_ca

    def CA_rmsd_residue_number(self, resi_vec):
        rms_ca = pyrosetta.rosetta.core.scoring.CA_rmsd(self.native_pose, self.pose, resi_vec)
        return rms_ca

    def all_atom_rmsd(self):
        rms_all = pyrosetta.rosetta.core.scoring.all_atom_rmsd(self.native_pose, self.pose)
        return rms_all
