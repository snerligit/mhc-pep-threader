#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: April 16, 2017
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *
from collections import defaultdict
# import rosetta files


# additional bio libraries

#custom libraries
from chain_split import CHAIN_SPLIT
from interface_analyzer import INTERFACE

# import other required libraries
import os
import sys
import argparse

if __name__ == "__main__":

    energies = defaultdict(dict)

    parser = argparse.ArgumentParser(description='Perform threading of the template structure onto the target sequence')
    parser.add_argument("-dir", help="provide dir name containing pdb")
    args = parser.parse_args()

    dirname = args.dir
    init()
    for files in os.listdir(dirname):
        if ".pdb" in files and "clean" not in files and "with" in files:
            my_pose = pose_from_pdb(dirname+"/"+files)
            chain = CHAIN_SPLIT(my_pose, 180)
            chain.cut()
            split = chain.get_pose()
            ia = INTERFACE(split)
            ia.analyze()
            fields = files.split("_")
            if len(fields) > 5:
                if fields[4] in energies:
                    energies[fields[4]] += ia.get_dG()
                else:
                    energies[fields[4]] = ia.get_dG()

    for key, value in energies.items():
        print(key+": "+str(value/100.0))
