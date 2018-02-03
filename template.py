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


class TEMPLATE:

    template_pdb = None
    templare_pose = None

    def __init__(self, pdb):
        self.template_pdb = pdb

    def get_pose(self):
        return self.template_pose

    def get_seq(self):
        return self.template_pose.sequence()

    def get_pdb(self):
        return self.template_pdb

    def get_name(self):
        template_pdb_fields = self.template_pdb.split(".")
        template_pdb_name = ".".join(template_pdb_fields[0:len(template_pdb_fields)-1])
        return template_pdb_name

    def get_stripped_name(self):
        template_pdb_name = os.path.basename(self.template_pdb)
        return template_pdb_name.split('.')[0]

    def get_template_path(self):
        return os.path.abspath(os.path.dirname(self.template_pdb))+"/"

    def append_clean_pdb(self):
        template_pdb_name = self.get_name()
        self.template_pdb = template_pdb_name + ".clean.pdb"

    def treat_template_structure(self, idealize_relax = False):
        if idealize_relax == True:
            cleanATOM(self.template_pdb)
            append_clean_pdb()
            self.template_pose = pose_from_pdb(template_clean_pdb)
            self.template_pose = IDEALIZE().idealize_pdb(template_pose)
            self.template_pose = RELAX().relax_pdb(template_pose)
        else:
            self.template_pose = pose_from_pdb(self.template_pdb)
