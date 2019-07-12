#! /home/snerli/anaconda3/bin/python
# additional bio libraries

from pyrosetta import *

#custom libraries
from scoring.rmsd import RMSD
from ia.chain_split import CHAIN_SPLIT
from ia.interface_analyzer import INTERFACE

# import other required libraries
import os
import sys
import fcntl
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

def rmsd(my_pose, args):
    native_pose = pose_from_pdb(args.native)

    rms_value = 0
    rms_obj = RMSD(my_pose, native_pose)

    rms_all = rms_obj.all_atom_rmsd()
    rms_ca = rms_obj.CA_rmsd()
    resi_list = extract_residues_from_movemap(args.movemap)
    rms_ca_movemap = rms_obj.CA_rmsd_residue_number(resi_list)

    return [rms_all, rms_ca, rms_ca_movemap]

def ia(my_pose, args):
    chain = CHAIN_SPLIT(my_pose, args.cutpoint)
    chain.cut()
    split = chain.get_pose()

    ia = INTERFACE(split)
    ia.analyze()
    return ia.get_dG()

def write_to_file(args, tag, bind_energy, rms_all, rms_ca, rms_ca_movemap):
    with open(args.outfile, "a") as writefilehandler:
        fcntl.flock(writefilehandler, fcntl.LOCK_EX)
        writefilehandler.write(tag+","+str(bind_energy)+","+str(rms_all)+","+str(rms_ca)+","+str(rms_ca_movemap)+"\n")
        fcntl.flock(writefilehandler, fcntl.LOCK_UN)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='compute rmsd values')
    parser.add_argument("-pdb", help="provide the pdb files")
    parser.add_argument("-native", help="provide native pdb file")
    parser.add_argument("-movemap", help="find rmsd value for the residues used in the movemap")
    parser.add_argument("-cutpoint", help="provide cutpoint", type=int)
    parser.add_argument("-outfile", help="provide output file name", default="rmsds_ia.csv")
    args = parser.parse_args()

    init()

    my_pose = pose_from_pdb(args.pdb)
    [rms_all, rms_ca, rms_ca_movemap] = rmsd(my_pose, args)
    bind_energy = ia(my_pose, args)
    tag = args.pdb.split(".pdb")[0]
    write_to_file(args, tag, bind_energy, rms_all, rms_ca, rms_ca_movemap)
