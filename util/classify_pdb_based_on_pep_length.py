#!~/anaconda3/bin/python

# import liraries for reading files in a directory
import os
import sys
import argparse

# library that allows messing with structures in PDB or mmcif format
from Bio.PDB import *

#def read_pdb():


# This is where the program starts
if __name__ == "__main__":

    # Parse the input arguments
    parser = argparse.ArgumentParser(description='Classify the MHC peptides into groups based on peptide length')
    parser.add_argument("-dir", help="provide the dir where the pdbs are present")
    args = parser.parse_args()

    dir_ = args.dir

    # loop over each pdb file
    for filename in os.listdir(dir_):
        new_filename = filename+".clean.pdb"
        os.system("grep ^ATOM "+filename+" > "+new_filename)
        if ".pdb" in new_filename:
            p = PDBParser()
            s = p.get_structure('X', dir_+"/"+new_filename)
            for model in s:
                for chain in model:
                    residue_count = 0
                    for residue in chain:
                        residue_count += 1
                    if residue.id[1] == residue_count:
                        if residue_count <= 15 and residue_count >= 7:
                            print (filename, residue_count)
                            out_dir = str(residue_count)+"mers"
                            if not os.path.isdir(dir_+"/"+out_dir):
                                os.system("mkdir "+out_dir)
                            os.system("mv "+filename+" "+out_dir)
                            break
