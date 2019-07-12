# additional bio libraries

from pyrosetta import *

#custom libraries
from rmsd import RMSD

# import other required libraries
import os
import sys
import argparse

def extract_residues_from_movemap(movemap):
    readfilehandler = open(movemap, "r")
    resi_list = pyrosetta.rosetta.std.list_unsigned_long_t()
    for line in readfilehandler:
        line = line.rstrip()
        if "RESIDUE" in line and "CHI" in line:
            resi = line.split(" ")[1]
            resi_list.append(int(resi))
    readfilehandler.close()
    return resi_list

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='compute rmsd values')
    parser.add_argument("-dir", help="provide directory containing pdb files")
    parser.add_argument("-native", help="provide native pdb file")
    parser.add_argument("-pattern", help="provide pattern to look for in the target file name")
    parser.add_argument("-movemap", help="find rmsd value for the residues used in the movemap")
    parser.add_argument("-choice", help="which type of rmsd do you want to calculate", choices=['all','CA','CA_with_start_and_end','CA_with_movemap'])
    parser.add_argument("-start", help="start residue index", type=int)
    parser.add_argument("-end", help="end residue index", type=int)
    parser.add_argument("-outfile", help="provide output file name", default="rmsds.txt")
    args = parser.parse_args()

    init()

    writefilehandler = open(args.outfile,"w")

    native_pose = pose_from_pdb(args.native)

    for filename in os.listdir(args.dir):
        if args.pattern in filename and ".pdb" in filename:
            rms_value = 0
            my_pose = pose_from_pdb(args.dir+"/"+filename)
            rms_obj = RMSD(my_pose, native_pose)

            if args.choice == "all":
                rms_value = rms_obj.all_atom_rmsd()
            elif args.choice == "CA":
                rms_value = rms_obj.CA_rmsd()
            elif args.choice == "CA_with_start_and_end":
                rms_value = rms_obj.CA_rmsd_between_indices(start, end)
            elif args.choice == "CA_with_movemap":
                resi_list = extract_residues_from_movemap(args.movemap)
                rms_value = rms_obj.CA_rmsd_residue_number(resi_list)
            writefilehandler.write(filename+" "+str(rms_value)+"\n")

    writefilehandler.close()
