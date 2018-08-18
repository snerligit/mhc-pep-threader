#!~/anaconda3/bin/python

# import liraries for reading files in a directory
import os
import sys
import argparse

from collections import defaultdict

# This is where the program starts
if __name__ == "__main__":

    # Parse the input arguments
    parser = argparse.ArgumentParser(description='Classify the MHC peptides into groups based on peptide length')
    parser.add_argument("-sc", help="provide the score file")
    parser.add_argument("-iedb", help="provide iedb file")
    args = parser.parse_args()

    sc_ = args.sc
    iedb_ = args.iedb

    ia_scores = defaultdict(list)
    iedb_scores = defaultdict(list)

    sc_file_handler = open(sc_, "r")
    for line in sc_file_handler:
        line = line.rstrip()
        fields = line.split()
        key = fields[1].split("_")[0]
        ia_scores[key].append(fields[0])
    sc_file_handler.close()

    iedb_file_handler = open(iedb_, "r")
    for line in iedb_file_handler:
        line = line.rstrip()
        fields = line.split(",")
        key = fields[0].split(".")[0]
        iedb_scores[key].append(line)
    iedb_file_handler.close()

    for key, value in iedb_scores.items():
        if key in ia_scores:
            print(iedb_scores[key][0]+","+ia_scores[key][0])
