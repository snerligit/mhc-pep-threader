#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: April 16 2, 2018
#   Email: snerli@ucsc.edu
#


import os
import sys
import subprocess
import argparse

aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def generate_peptides(peptide_sequence, type, skip_N=1, skip_C=8):

    if type == 'ala_scan':
        amino_acids = list(peptide_sequence)
        for aa in range(0,len(amino_acids)):
            amino_acids[aa] = 'A'
            print (">"+"".join(amino_acids))
            print ("".join(amino_acids))
            #print ("human   HLA A*0201              9       1       "+"".join(amino_acids)+"       >       20000.0")
            amino_acids = list(peptide_sequence)
    elif type == 'pos_saturated':
        amino_acids = list(peptide_sequence)
        for aa in range(0,len(amino_acids)):
            acutal_aa = amino_acids[aa]
            if aa == skip_N or aa == skip_C:
                for a in aa_list:
                    amino_acids[aa] = a
                    if (acutal_aa != a):
                        print (">"+"".join(amino_acids))
                        print ("".join(amino_acids))
                        amino_acids = list(peptide_sequence)
    else:
        print("Invalid type choice. Available options are \"ala_scan\" and \"pos_saturated\"")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create peptide sequnce library')
    parser.add_argument("-peptide_sequence", help="provide the peptide sequence", required=True)
    parser.add_argument("-type", choices=['ala_scan','pos_saturated'], default='ala_scan', help="type of peptides to generate")


    args = parser.parse_args()

    peptide_sequence_ = args.peptide_sequence
    type_ = args.type

    generate_peptides(peptide_sequence_, type_)
