# additional bio libraries

from pyrosetta import *

#custom libraries
from chain_split import CHAIN_SPLIT
from interface_analyzer import INTERFACE

# import other required libraries
import os
import sys
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='call ia')
    parser.add_argument("-dir", help="provide directory containing pdb names")
    parser.add_argument("-cutpoint", help="provide cutpoint", type=int)
    parser.add_argument("-pattern", help="provide pattern to look for in the file name", default="relaxed")
    parser.add_argument("-outfile", help="provide output file name")
    args = parser.parse_args()

    init()

    writefilehandler = open(args.outfile,"w")

    for filename in os.listdir(args.dir):
        if args.pattern in filename and ".pdb" in filename:
            my_pose = pose_from_pdb(args.dir+"/"+filename)
            chain = CHAIN_SPLIT(my_pose, args.cutpoint)
            chain.cut()
            split = chain.get_pose()

            ia = INTERFACE(split)
            ia.analyze()
            print("Interface energy: "+str(ia.get_dG()))
            writefilehandler.write(filename+" "+str(ia.get_dG()))

    writefilehandler.close()
