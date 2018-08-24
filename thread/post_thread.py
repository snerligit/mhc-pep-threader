#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: May 25, 2017
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

# additional bio libraries

#custom libraries
from idealize_relax.relax import RELAX
from idealize_relax.movemap import MOVEMAP
from ia.chain_split import CHAIN_SPLIT
from ia.interface_analyzer import INTERFACE

# import other required libraries
import os
import sys

class POST_THREADING:

    args = None
    threaded_pose = None
    tag = ""
    pre_threader = None
    movemap = None
    mpi_install = True

    def __init__(self, pre_threader, threaded_pose, tag):

        self.args = pre_threader.args
        self.pre_threader = pre_threader
        self.threaded_pose = threaded_pose
        self.tag = tag

    def treatment_post_homology_modeling(self):

        self.threaded_pose.dump_pdb(self.tag+".pdb")
        self.movemap = MOVEMAP(self.tag+".pdb", self.args.get_pep_start_index(),
                            self.pre_threader.get_pep_length(), self.args.get_groove_distance(),
                            self.tag+".movemap", True)
        self.movemap.apply()

        if self.args.get_relax_after_threading():

            relaxed_pose = self.threaded_pose

            # yet to make this work
            try:
                from jd.job_distributor import JOB_DISTRIBUTOR
                job_dist = JOB_DISTRIBUTOR()
                print ("No. of jobs for relax: ", self.args.get_nstruct())
                job_dist.SimpleMPIJobDistributor(self.args.get_nstruct(), self.minimize_and_calculate_energy)
            except ImportError:
                self.mpi_install = False

            #self.mpi_install = False
            if not self.mpi_install:
                for job_id in range(self.args.get_nstruct()):
                    self.minimize_and_calculate_energy(job_id)

    def minimize_and_calculate_energy(self, i):

        relax = RELAX()
        threaded_pose = Pose()
        threaded_pose.detached_copy(self.threaded_pose)

        relaxed_threaded_pose = Pose()
        relaxed_threaded_pose = relax.relax_pdb_with_movemap(threaded_pose, self.movemap.get_movemap())
        relaxed_threaded_pose.dump_pdb(self.tag+"_relaxed_"+str(i)+".pdb")

        chain_split = CHAIN_SPLIT(relaxed_threaded_pose, self.args.get_interface_cupoint())
        chain_split.cut()
        split_pose = Pose()
        split_pose = chain_split.get_pose()

        ia = INTERFACE(split_pose)
        ia.analyze()
        print("Interface energy: ", ia.get_dG())
