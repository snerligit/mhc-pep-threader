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
from collections import defaultdict


aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def generate_peptides(peptide_sequence, type, skip_anchors, skip_N=1, skip_C=9):

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

        print (">"+peptide_sequence)
        print (peptide_sequence)

        for aa in range(0,len(amino_acids)):
            amino_acids[aa] = 'A'

        print (">"+"".join(amino_acids))
        print ("".join(amino_acids))

        amino_acids = list(peptide_sequence)
        for aa in range(0,len(amino_acids)):
            acutal_aa = amino_acids[aa]
            if skip_anchors == True:
                if aa != skip_N and aa != skip_C:
                    for a in aa_list:
                        amino_acids[aa] = a
                        if (acutal_aa != a):
                            print (">"+"".join(amino_acids))
                            print ("".join(amino_acids))
                            amino_acids = list(peptide_sequence)
            else:
                for a in aa_list:
                    amino_acids[aa] = a
                    if (acutal_aa != a):
                        print (">"+"".join(amino_acids))
                        print ("".join(amino_acids))
                        amino_acids = list(peptide_sequence)
    elif type == 'pos_saturated_double':
        amino_acids = list(peptide_sequence)

        peptide_sequences = defaultdict(dict)

        peptide_sequences[peptide_sequence] = 1

        for aa in range(0,len(amino_acids)):
            amino_acids[aa] = 'A'

        peptide_sequences["".join(amino_acids)] = 1

        amino_acids = list(peptide_sequence)
        ctr = 0
        for aa1 in range(0,len(amino_acids)):
            actual_aa1 = amino_acids[aa1]
            aa1_index = aa1
            if aa1_index != skip_N and aa1_index != skip_C:
                for a1 in aa_list:
                    if actual_aa1 != a1:
                        amino_acids[aa1_index] = a1
                        for aa2 in range(0,len(amino_acids)):
                            aa2_index = aa2
                            actual_aa2 = amino_acids[aa2]
                            if aa2_index != skip_C and aa2_index != skip_N and aa2 != aa1:
                                for a2 in aa_list:
                                    if actual_aa2 != a2:
                                        amino_acids[aa2_index] = a2
                                        my_seq = "".join(amino_acids)
                                        peptide_sequences[my_seq] = 1
                                        amino_acids[aa2] = actual_aa2
            amino_acids = list(peptide_sequence)


        for key, value in peptide_sequences.items():
            print ">"+key
            print key
            a = 10

    else:
        print("Invalid type choice. Available options are \"ala_scan\" and \"pos_saturated\"")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create peptide sequnce library')
    parser.add_argument("-peptide_sequence", help="provide the peptide sequence", required=True)
    parser.add_argument("-type", choices=['ala_scan','pos_saturated', 'pos_saturated_double'], default='ala_scan', help="type of peptides to generate")
    parser.add_argument("-skip_anchors", help="set to true if you don't want to change the anchors",  action='store_true')

    args = parser.parse_args()

    peptide_sequence_ = args.peptide_sequence
    type_ = args.type
    skip_anchors_ = args.skip_anchors

    generate_peptides(peptide_sequence_, type_, skip_anchors_)
