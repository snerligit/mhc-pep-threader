#!/usr/bin/python

#       Sgourakis Lab
#   Author: Hailey Wallace
#   Date: June 21, 2019
#   Email: hmwalla2@ncsu.edu
#   Edits to allow for backbone perturbations before FastRelax
#

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

#custom libraries
from idealize_relax.movemap import MOVEMAP
from ia.chain_split import CHAIN_SPLIT
from ia.interface_analyzer import INTERFACE
from input_output.output.output_interface_energies import OUT_INTERFACE_ENERGY

from idealize_relax.backrub import BACKRUB
from idealize_relax.relax import RELAX

# import other required libraries
import os
import sys

'''

POST_THREADING class contains all the necessary functionalities required to
perform operations post threading process.

'''

class POST_THREADING:

    # class members
    args = None # store the input arguments
    threaded_pose = None # store the threaded pose
    tag = "" # tag is used for naming the output file
    pre_threader = None # pre threader object
    movemap = None # movemap object
    mpi_install = True # should I run mpi?
    output_energies = None # output the interface energies

    # constructor
    def __init__(self, pre_threader, threaded_pose, tag):
        self.args = pre_threader.args
        self.pre_threader = pre_threader
        self.threaded_pose = threaded_pose
        self.tag = tag
        self.output_energies = OUT_INTERFACE_ENERGY(self.args.get_out_file())

    # method to refine the structures post homology modeling
    def treatment_post_homology_modeling(self):
        # output threaded structure
        self.threaded_pose.dump_pdb(self.tag+".pdb")
        # create movemap
        self.movemap = MOVEMAP(self.tag+".pdb", self.args.get_pep_start_index(),
                            self.pre_threader.get_pep_length(), self.args.get_groove_distance(),
                            self.tag+".movemap")
        self.movemap.apply()

        # if refine after relax
        if self.args.get_relax_after_threading():

            relaxed_pose = self.threaded_pose

            # Can we split this into multiple nodes?
            # if we find more
            try:
                from jd.job_distributor import JOB_DISTRIBUTOR
                job_dist = JOB_DISTRIBUTOR()
                print ("No. of jobs for relax: ", self.args.get_nstruct())
                job_dist.SimpleMPIJobDistributor(self.args.get_nstruct(), self.minimize_and_calculate_energy)
            except ImportError:
                self.mpi_install = False

            # if mpi is not installed, just run on a
            # single thread
            if not self.mpi_install:
                for job_id in range(self.args.get_nstruct()):
                    self.minimize_and_calculate_energy(job_id)
            self.output_energies.write()

    #method to perturb backbone
    def backbone_mover(self, i):
        backrub = BACKRUB()
        threaded_pose = Pose()
        threaded_pose.detached_copy(self.threaded_pose)


    # method to perform energy minimization and calculate
    # interface energies
    def minimize_and_calculate_energy(self, i):
        relax = RELAX()
        threaded_pose = Pose()
        threaded_pose.detached_copy(self.threaded_pose)


        # output minimized structures
        relaxed_threaded_pose = Pose()
        relaxed_threaded_pose = relax.relax_pdb_with_movemap(threaded_pose, self.movemap.get_movemap())
        relaxed_threaded_pose.dump_pdb(self.tag+"_relaxed_"+str(i)+".pdb")

        freshly_relaxed_pose = pose_from_pdb(self.tag+"_relaxed_"+str(i)+".pdb")
        # split the chains of the relaxed structures
        # according to the cutpoint because threader
        # create a singe chain
        chain_split = CHAIN_SPLIT(freshly_relaxed_pose, self.args.get_interface_cupoint())
        chain_split.cut()
        split_pose = Pose()
        split_pose = chain_split.get_pose()

        # calculate interface energies
        ia = INTERFACE(split_pose)
        ia.analyze()
        print("Interface energy: ", ia.get_dG())
        self.output_energies.add(self.tag, ia.get_dG())
