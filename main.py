#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: February 2, 2018
#   Email: snerli@ucsc.edu
#


# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

#custom libraries
from model_hla import MODEL_HLA
from template import TEMPLATE
from database.HLA_sequences_180 import hla_sequences

# import other required libraries
import os
import sys
import subprocess
import argparse

def model_hla_for_each_peptide(template_pdb_, template_, fasta_, groove_distance_, nstruct_, movemap_, peptide_, mhcs_, beta2m):

    if beta2m == None:
        beta2m = ""
    pep_file_handler = open(peptide_, "r")
    for line in pep_file_handler:
        line = line.rstrip()
        if ">" in line:
            pep_header = line[1:]
        else:
            pep_seq = line
            modeller = MODEL_HLA(template_pdb_, template_, fasta_, groove_distance_, nstruct_, movemap_, pep_header, pep_seq, mhcs_, beta2m)
            modeller.apply()
    pep_file_handler.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Perform threading of the template structure onto the target sequence')
    parser.add_argument("-template_pdb", help="provide template structure in PDB to perform threading")
    parser.add_argument("-fasta", help="provide FASTA format file of the mhc alleles to be modeled")
    parser.add_argument("-peptide", help="provide fasta file with peptide sequences that need to be threaded")
    parser.add_argument("-movemap", help="provide a movemap file for relaxing the target")
    parser.add_argument("-groove_distance", help="provide distance to select nearest groove residues from the peptide", type=float, default=3.5)
    parser.add_argument("-idealize_relax", help="idealize and relax template structure before threading", action='store_true')
    parser.add_argument("-nstruct", help="number of times a threaded structure should be relaxed", type=int, default=1)
    parser.add_argument("-mhcs", help="provide the list of names of MHCs in the file, if you want to include all, just type \'all\' in the file")
    parser.add_argument("-list_mhcs", help="List all the HLAs for which sequences are available in the database", action='store_true')
    parser.add_argument("-beta2m", help="provide beta2m sequence")
    parser.add_argument("-mhc_chain", help="provide mhc chain id in the template")
    parser.add_argument("-peptide_chain", help="provide peptide chain id in the template")

    args = parser.parse_args()

    template_pdb_ = args.template_pdb
    fasta_ = args.fasta
    idealize_relax_ = args.idealize_relax
    nstruct_ = args.nstruct
    movemap_ = args.movemap
    groove_distance_ = args.groove_distance
    peptide_ = args.peptide
    list_mhcs_ = args.list_mhcs
    mhcs_ = args.mhcs
    beta2m_ = args.beta2m
    mhc_chain_ = args.mhc_chain
    peptide_chain_ = args.peptide_chain

    if list_mhcs_ == True:
        for key,value in hla_sequences.items():
            print(key)
    else:
        if (fasta_ == None and mhcs_ == None) or template_pdb_ == None or peptide_ == None:
            print("Please provide mhc list, peptide list and template_pdb to perform threading\n Type -h to see all the options")
            exit(1)
        # Load Rosetta datbase files
        init()
        template_ = TEMPLATE(template_pdb_)
        template_.treat_template_structure(mhc_chain_, peptide_chain_)
        model_hla_for_each_peptide(template_pdb_, template_, fasta_, groove_distance_, nstruct_, movemap_, peptide_, mhcs_, beta2m_)
