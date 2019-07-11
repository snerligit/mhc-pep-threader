
from pyrosetta import *

import os
import sys
import argparse

from movemap import MOVEMAP

def exists(tag, dir):

    for filename in os.listdir(dir):
        if ".movemap" in filename:
            if tag == filename.split(".movemap")[0]:
                return True

    return False

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='generate movemap')
    parser.add_argument("-pdb", help="provide pdb name", default=None)
    parser.add_argument("-dir", help="provide directory containing pdbs", default=None)
    parser.add_argument("-pattern", help="provide a pattern", default="low")
    parser.add_argument("-pep_start_index", help="provide peptide start index", type=int)
    parser.add_argument("-pep_length", help="provide peptide length", type=int)
    parser.add_argument("-groove_distance", help="provide groove distance", type=float, default=3.5)
    args = parser.parse_args()

    init()

    if (args.pdb == None and args.dir != None):
        for filename in os.listdir(args.dir):
            if args.pattern in filename and ".movemap" not in filename:
                tag = args.dir+"/"+filename.split(".pdb")[0]
                if not exists(tag, args.dir):
                    movemap = MOVEMAP(tag+".pdb", args.pep_start_index, args.pep_length, args.groove_distance, tag+".movemap")
                    movemap.apply()
                else:
                    print ("Ignored: ", tag)
    elif args.pdb != None:
        tag = args.pdb.split(".pdb")[0]
        movemap = MOVEMAP(tag+".pdb", args.pep_start_index, args.pep_length, args.groove_distance, tag+".movemap")
        movemap.apply()
