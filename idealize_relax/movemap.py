#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: December 1, 2017
#   Email: snerli@ucsc.edu
#

'''

MOVEMAP class contains all the necessary functionalities required to create a
a set of residues whose backbone and side-chain torsion angles are allowed to
modify during structure refinement.

'''

# import other required libraries
import sys
import math

# import biopython libraries
#
# Hamelryck, T., Manderick, B. (2003) PDB parser and structure class implemented in Python. Bioinformatics 19: 2308-2310
from Bio.PDB import *
#

# Rosetta libraries
from pyrosetta.rosetta.core.kinematics import MoveMap

# custom libraries
from util.custom_pdb import MyPDB


class MOVEMAP:

    # class members
    structure = None # store the pdb file
    pep_start_index = None # store the starting residue index from which you want to create movemap
    pep_len = -1 # length of a peptide
    groove_distance = 3.5 # also include those residues that are at a distance of 3.5 A from peptide
    filename = "" # movemap file name

    protein_residues = []
    peptide_residues = []
    movemap_res_list = []

    # constructor
    def __init__(self, pdb_file, pep_start_index, pep_length, groove_distance, filename):
        self.structure = MyPDB(pdb_file)
        self.structure.read_chain_A() # This should work since movemap is created after threading is done which creates only a single chain A
        self.pep_start_index = pep_start_index
        self.groove_distance = groove_distance
        self.filename = filename
        self.pep_len = pep_length

    # method to create a write to movemap file
    def apply(self):
        self.protein_residues = []
        self.peptide_residues = []
        self.movemap_res_list = []
        self.extract_peptide_protein_residues()
        self.create_pep_protein_list()
        self.write_movemap()

    # method to make a list of all peptide and protein residues that
    # will be filtered in the next stage
    def extract_peptide_protein_residues(self):
        residue_number = 1
        print (self.pep_start_index+self.pep_len)
        for residue in self.structure.get_pdb():
            if residue_number >= self.pep_start_index and residue_number < self.pep_start_index+self.pep_len:
                self.peptide_residues.append(residue)
            else:
                self.protein_residues.append(residue)
            residue_number += 1

    # method to fetch all the residues that can be moved during refinement
    def create_pep_protein_list(self):
        for pep_res in self.peptide_residues: # for each peptide residue
            for protein_res in self.protein_residues: # for each protein residue
                for pep_atom in pep_res: # for each atom in the pep residue
                    for protein_atom in protein_res: # for each atom in the protein residue
                        distance = math.fabs(pep_res[pep_atom.get_name()] - protein_res[protein_atom.get_name()]) # calculate distance between two atoms
                        if distance <= self.groove_distance:
                            residue_number = int(str(protein_res).split()[3].split("=")[1])
                            if residue_number not in self.movemap_res_list:
                                self.movemap_res_list.append(residue_number)
            self.movemap_res_list.append(int(str(pep_res).split()[3].split("=")[1]))

    # method to write to a movemap file
    # see here for movemap file requirements: https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/movemap-file
    def write_movemap(self):
        movemap_handler = open(self.filename, "w")
        movemap_handler.write("RESIDUE * NO\n")
        for residue in self.movemap_res_list:
            movemap_handler.write("RESIDUE "+str(residue)+" BBCHI\n")
        movemap_handler.write("JUMP * YES\n")
        movemap_handler.close()

    # getter methods
    def get_filename(self):
        return self.filename

    # method to create Rosetta specific movemap object
    # read from file and return thr object
    def get_movemap(self):
        movemap_object = MoveMap()
        movemap_object.init_from_file(self.filename)
        return movemap_object

    def extract_residues_from_movemap(movemap):
        readfilehandler = open(movemap, "r")
        resi_list = pyrosetta.rosetta.std.list_unsigned_long_t()
        for line in readfilehandler:
            line = line.rstrip()
            if "RESIDUE" in line and "CHI" in line:
                resi = line.split(" ")[1]
                resi_list.append(int(resi))
        readfilehandler.close()
        return resi_list
