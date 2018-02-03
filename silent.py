#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: February 2, 2018
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *
from pyrosetta.rosetta.core.io.silent import *

class SILENT:

    silent_file_data = None
    silent_struct = None
    silent_struct_type = "binary"

    def __init__(self):
        sf_options = SilentFileOptions()
        self.silent_file_data = SilentFileData(sf_options)
        self.silent_struct = SilentStructFactory.get_instance().get_silent_struct(self.silent_struct_type)

    def add(self, pose, tag, filename):
        self.silent_struct.fill_struct(pose, tag)
        self.silent_file_data.write_silent_struct(self.silent_struct, filename)
