#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: December 1, 2017
#   Email: snerli@ucsc.edu
#

# import other required libraries
import os
import sys
import math
import numpy
import argparse
import subprocess

# import biopython libraries
#
# Hamelryck, T., Manderick, B. (2003) PDB parser and structure class implemented in Python. Bioinformatics 19: 2308-2310
from Bio.PDB import *
#

# Rosetta libraries
from pyrosetta.rosetta.core.kinematics import MoveMap

# custom libraries
from custom_pdb import MyPDB


class MOVEMAP:

    structure = None
    pep_start_index = None
    groove_distance = 3.5
    filename = ""
    already_computed = False

    protein_residues = []
    peptide_residues = []
    movemap_res_list = []

    def __init__(self, pdb_file, pep_start_index, groove_distance, filename, already_computed=False):
        self.structure = MyPDB(pdb_file)
        self.structure.read_chain_A()
        self.pep_start_index = pep_start_index
        self.groove_distance = groove_distance
        self.filename = filename
        self.already_computed = already_computed

    def apply(self):
        self.extract_peptide_protein_residues()
        self.create_pep_protein_list()
        self.write_movemap()

    def get_movemap(self):
        if self.already_computed == False:
            self.apply()
        movemap_object = MoveMap()
        movemap_object.init_from_file(self.filename)
        return movemap_object

    def write_movemap(self):
        movemap_handler = open(self.filename, "w")
        movemap_handler.write("RESIDUE * NO\n")
        for residue in self.movemap_res_list:
            movemap_handler.write("RESIDUE "+str(residue)+" BBCHI\n")
        movemap_handler.write("JUMP * YES\n")
        movemap_handler.close()

    def create_pep_protein_list(self):
        for pep_res in self.peptide_residues:
            for protein_res in self.protein_residues:
                for pep_atom in pep_res:
                    for protein_atom in protein_res:
                        distance = math.fabs(pep_res[pep_atom.get_name()] - protein_res[protein_atom.get_name()])
                        if distance <= self.groove_distance:
                            residue_number = int(str(protein_res).split()[3].split("=")[1])
                            if residue_number not in self.movemap_res_list:
                                self.movemap_res_list.append(residue_number)
            self.movemap_res_list.append(int(str(pep_res).split()[3].split("=")[1]))

    def extract_peptide_protein_residues(self):
        residue_number = 1
        for residue in self.structure.get_pdb():
            if residue_number >= self.pep_start_index:
                self.peptide_residues.append(residue)
            else:
                self.protein_residues.append(residue)
            residue_number += 1
