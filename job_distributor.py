#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 31, 2018
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *
from pyrosetta.rosetta.core.scoring import *

class JOB_DISTRIBUTOR:

    filename = ""
    nstruct = 1
    score_function = None
    job_dist_obj = None

    def __init__(self, filename, nstruct=1):
        self.filename = filename
        #self.score_function = score_function
        self.score_function = get_fa_scorefxn()
        self.nstruct = nstruct
        self.job_dist_obj = PyJobDistributor(self.filename, self.nstruct, self.score_function)

    def get_dist_obj(self):
        return self.job_dist_obj
