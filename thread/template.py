#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 26, 2017
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *
from pyrosetta.toolbox import cleanATOM

# import rosetta files
from pyrosetta.rosetta.protocols.relax import *
from pyrosetta.rosetta.core.scoring import *

from Bio.PDB import *

from idealize_relax.idealize import IDEALIZE
from idealize_relax.relax import RELAX


class TEMPLATE:

    template_pdb = None
    template_pose = None

    def __init__(self, pdb):
        self.template_pdb = pdb

    def get_pose(self):
        return self.template_pose

    def get_sequence(self):
        return self.template_pose.sequence()

    def get_pdb(self):
        return self.template_pdb

    def get_pose_from_pdb(self):
        if self.template_pose == None:
            self.template_pose = pose_from_pdb(self.template_pdb)
        return self.template_pose

    def get_name(self):
        template_pdb_fields = self.template_pdb.split(".")
        template_pdb_name = ".".join(template_pdb_fields[0:len(template_pdb_fields)-1])
        return template_pdb_name

    def get_stripped_name(self):
        template_pdb_name = os.path.basename(self.template_pdb)
        return template_pdb_name.split('.')[0]

    def get_template_path(self):
        return os.path.abspath(os.path.dirname(self.template_pdb))+"/"

    def trim_pdb(self, mhc_chain, peptide_chain, mhc_trim_length):
        p = PDBParser()
        s = p.get_structure('X', self.template_pdb)
        class ResSelect(Select):
            def accept_residue(self, res):
                if res.id[1] >= mhc_trim_length and res.parent.id == mhc_chain:
                    return False
                if res.parent.id != mhc_chain and res.parent.id != peptide_chain:
                    return False
                else:
                    return True
        io = PDBIO()
        io.set_structure(s)
        self.append_cropped_pdb()
        io.save(self.template_pdb, ResSelect())

    def save_pdb(self):
        self.template_pose.dump_pdb(self.template_pdb)

    def append_clean_pdb(self):
        template_pdb_name = self.get_name()
        self.template_pdb = template_pdb_name + ".clean.pdb"

    def append_cropped_pdb(self):
        template_pdb_name = self.get_name()
        self.template_pdb = template_pdb_name + ".cropped.pdb"

    def treat_template_structure(self, mhc_chain, peptide_chain, trim_mhc, mhc_trim_length, idealize_relax = False):

        cleanATOM(self.template_pdb)
        self.append_clean_pdb()
        if trim_mhc:
            self.trim_pdb(mhc_chain, peptide_chain, mhc_trim_length)
        self.template_pose = pose_from_pdb(self.template_pdb)
        if idealize_relax == True:
            template_pose = self.template_pose
            self.template_pose = IDEALIZE().idealize_pdb(template_pose)
            self.template_pose = RELAX().relax_pdb(template_pose)
        self.save_pdb()
