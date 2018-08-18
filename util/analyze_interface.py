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
        #new_filename = filename+".clean.pdb"
        #os.system("grep ^ATOM "+filename+" > "+new_filename)
        new_filename = filename
        if ".pdb" in new_filename and "relaxed" in new_filename:
            p = PDBParser()
            s = p.get_structure('X', dir_+"/"+new_filename)
            pep_chain = "C"
            for model in s:
                for chain in model:
                    print(chain.id)
                    if chain.id == "P":
                        pep_chain = "P"
            os.system("/home/snerli/rosetta/Rosetta_ref2015/main/source/bin/InterfaceAnalyzer.linuxgccrelease -database /home/snerli/rosetta/Rosetta_ref2015/main/database -s "+new_filename+" -out:file:score_only "+new_filename+"_ia.sc -interface AB_"+pep_chain)
