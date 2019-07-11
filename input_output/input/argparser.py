#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: May 25, 2018
#   Email: snerli@ucsc.edu
#

'''

ARGPARSE class contains all the necessary functionalities required to accept input from
users.

'''

import argparse


class ARGPARSE:

    # class member
    args = None # stores all the arguments

    # constructor
    def __init__(self):
        self.parse_args()
        self.check_user_input()

    # method that provides parser arguments that are seen on the command line when
    # running this method
    def parse_args(self):
        parser = argparse.ArgumentParser(description='Perform threading of the template structure onto the target sequence')
        parser.add_argument("-list_mhcs", help="List all the HLAs for which sequences are available in the database", action='store_true')
        parser.add_argument("-template_pdb", help="provide template structure in PDB to perform threading")
        parser.add_argument("-mhcs", help="provide the list of names of MHCs in the file, if you want to include all, just type \'all\' in the file")
        parser.add_argument("-beta2m", help="provide the file containing names of beta2m, Choices include: choices=[humanbeta2m, mousebeta2m, chickenbeta2m, bovinebeta2m, ratbeta2m, none]")
        parser.add_argument("-peptides", help="provide fasta file with peptide sequences that need to be threaded")
        parser.add_argument("-tcr", help="provide the file containing tcr sequence in fasta format")
        parser.add_argument("-chaperone", help="provide the file containing names of chaperones, Choices include: choices=[tapasin, tapbpr, none]")
        parser.add_argument("-mhc_trim_length", help="provide the number of residue from which the mhc should be trimmed", type=int, default=181)
        parser.add_argument("-no_trim_mhc", help="Should we model the whole complex", action='store_false')
        parser.add_argument("-idealize_relax", help="idealize and relax template structure before threading", action='store_true')
        parser.add_argument("-relax_after_threading", help="idealize and relax template structure before threading", action='store_true')
        parser.add_argument("-mhc_chain", help="provide mhc chain id in the template")
        parser.add_argument("-peptide_chain", help="provide peptide chain id in the template")
        parser.add_argument("-pep_start_index", help="provide peptide start index", type=int, default=0)
        parser.add_argument("-groove_distance", help="provide distance to select nearest groove residues from the peptide", type=float, default=3.5)
        parser.add_argument("-interface_cutpoint", help="last residue index that separates the interfaces for which you are calculating binding energies", type=int, default=0)
        parser.add_argument("-out_file", help="output file name in csv format to write the binding energies", default="binding_energies.csv")
        parser.add_argument("-nstruct", help="number of times a threaded structure should be relaxed", type=int, default=1)
        parser.add_argument("-clustal_path", help="Path to clustal omega", default='clustalo')
        parser.add_argument("-native", help="native pdb file to compare and report RMSD values", default=None)

        self.args = parser.parse_args()

    # method to validate user specified input
    # check if a more extensive check is required
    def check_user_input(self):
        if self.is_list_mhcs() == False:
            if self.get_mhcs() == None or self.get_template_pdb() == None:
                print("Please provide mhc list, peptide list and template_pdb to perform threading\n Type -h to see all the options")
                exit(1)

    # methods that store actions
    def is_list_mhcs(self):
        return self.args.list_mhcs

    def is_no_trim_mhc_flag_set(self):
        return self.args.no_trim_mhc

    # getter methods that return all user provided inputs
    def get_template_pdb(self):
        return self.args.template_pdb

    def get_interface_cupoint(self):
        return self.args.interface_cutpoint

    def get_mhc_trim_length(self):
        return self.args.mhc_trim_length

    def get_mhcs(self):
        return self.args.mhcs

    def get_peptides(self):
        return self.args.peptides

    def get_groove_distance(self):
        return self.args.groove_distance

    def get_nstruct(self):
        return self.args.nstruct

    def get_idealize_relax(self):
        return self.args.idealize_relax

    def get_relax_after_threading(self):
        return self.args.relax_after_threading

    def get_beta2m(self):
        if self.is_no_trim_mhc_flag_set():
            return None
        return self.args.beta2m

    def get_tcr(self):
        return self.args.tcr

    def get_chaperone(self):
        if self.is_no_trim_mhc_flag_set():
            return None
        return self.args.chaperone

    def get_peptide_chain(self):
        return self.args.peptide_chain

    def get_mhc_chain(self):
        return self.args.mhc_chain

    def get_pep_start_index(self):
        return self.args.pep_start_index

    def get_out_file(self):
        return self.args.out_file

    def get_clustal_path(self):
        return self.args.clustal_path

    def get_native(self):
        return self.args.native
