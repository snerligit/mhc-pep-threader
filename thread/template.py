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
from idealize_relax.idealize import IDEALIZE
from idealize_relax.relax import RELAX

# additional bio library
from Bio.PDB import *

'''

TEMPLATE class contains all the necessary functionalities required to
perform prepare the template structure for threading.

'''


class TEMPLATE:

    # class members
    template_pdb = None # store template pdb file name
    template_pose = None # store the template pose

    # constructor
    def __init__(self, pdb):
        self.template_pdb = pdb

    # method to return pose
    def get_pose(self):
        return self.template_pose

    # method to return sequence of the template
    def get_sequence(self):
        return self.template_pose.sequence()

    # method to return name of the template pdb file
    def get_pdb(self):
        return self.template_pdb

    # method to load pose from a pdb file
    def get_pose_from_pdb(self):
        if self.template_pose == None:
            self.template_pose = pose_from_pdb(self.template_pdb)
        return self.template_pose

    # method to return pdb name fields joined with dot
    # some of the naming convention we use
    def get_name(self):
        template_pdb_fields = self.template_pdb.split(".")
        template_pdb_name = ".".join(template_pdb_fields[0:len(template_pdb_fields)-1])
        return template_pdb_name

    # method to get the striped pdb name or file name without the .pdb extension
    def get_stripped_name(self):
        template_pdb_name = os.path.basename(self.template_pdb)
        return template_pdb_name.split('.')[0]

    # methot to return absolute path of the file where template is present
    def get_template_path(self):
        return os.path.abspath(os.path.dirname(self.template_pdb))+"/"

    # method to trim the pdb if specified
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

    # method to save pdb file
    def save_pdb(self):
        self.template_pose.dump_pdb(self.template_pdb)

    # method to append cleaned pdb after cleaning the template
    def append_clean_pdb(self):
        template_pdb_name = self.get_name()
        self.template_pdb = template_pdb_name + ".clean.pdb"

    # method to append cropped pdb after cropping the template structure
    def append_cropped_pdb(self):
        template_pdb_name = self.get_name()
        self.template_pdb = template_pdb_name + ".cropped.pdb"

    # method to clean, crop, idealize or relax pdb file
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
