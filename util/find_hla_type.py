#!~/anaconda3/bin/python

# import liraries for reading files in a directory
import os
import re
import sys
import argparse

# library that allows messing with structures in PDB or mmcif format
from Bio.PDB import *

threeToOne= {'ALA':'A',
            'GLY':'G',
            'VAL':'V',
            'ILE':'I',
            'LEU':'L',
            'TRP':'W',
            'PHE':'F',
            'TYR':'Y',
            'THR':'T',
            'SER':'S',
            'MET':'M',
            'CYS':'C',
            'PRO':'P',
            'LYS':'K',
            'ARG':'R',
            'HIS':'H',
            'ASP':'D',
            'ASN':'N',
            'GLU':'E',
            'GLN':'Q'
}

def find_peptide_seq(filename):
    new_filename = filename+".clean.pdb"
    os.system("grep ^ATOM "+filename+" > "+new_filename)
    if ".pdb" in new_filename:
        p = PDBParser()
        s = p.get_structure('X', dir_+"/"+new_filename)
        for model in s:
            for chain in model:
                residue_count = 0
                residues = []
                for residue in chain:
                    residue_count += 1
                    residues.append(threeToOne[residue.get_resname()])
                if residue.id[1] == residue_count:
                    if residue_count <= 15 and residue_count >= 7:
                        return "".join(residues)

def get_ic50(mhc, pep):
    filehandler = open("test.pep", "w")
    filehandler.write(pep)
    filehandler.close()
    os.system("netMHCpan -p test.pep -a "+mhc+" -BA -xls -l "+str(len(pep)))
    os.system("cat NetMHCpan_out.xls | awk '{print $7}' | grep -v \"nM\" > ic50")
    filehandler = open("ic50", "r")
    ic50_value = -1
    for line in filehandler:
        line = line.rstrip()
        if len(list(line)) > 0:
            ic50_value = line
    filehandler.close()
    return ic50_value

# This is where the program starts
if __name__ == "__main__":

    # Parse the input arguments
    parser = argparse.ArgumentParser(description='Classify the MHC peptides into groups based on peptide length')
    parser.add_argument("-dir", help="provide the dir where the pdbs are present")
    args = parser.parse_args()

    dir_ = args.dir

    # loop over each pdb file
    result = []
    for filename in os.listdir(dir_):
        if ".pdb" in filename:
            filehandler = open(dir_+"/"+filename, "r")
            for line in filehandler:
                if "ATOM" in line:
                    break
                line = line.rstrip()
                if "HLA-A" in line:
                    fields = line.split()
                    key = ""
                    for f in fields:
                        if "HLA-A" in f:
                            key = f
                            new_key = re.sub('\*', '', key)
                            org_key = new_key
                            type = list(new_key.split("-")[1])
                            new_key = "HLA-"
                            new_type = "".join(type)
                            if "A2.1" in new_type:
                                new_type = "A02:01"
                            if len(type) == 1 and type[0] == "A":
                                new_type = "A02:01"
                            if len(type) == 2 and (type[1] == "2" or type[1] == "," or type[1] == ":" or type[1] == ";"):
                                new_type = "A02:01"
                            if len(type) == 2 and (type[1] == "1"):
                                new_type = "A01:01"
                            if len(type) == 3 and (type[1] == "2" or type[1] == "3") and (type[2] == ";" or type[2] == ":" or type[2] == ","):
                                new_type = "A0"+type[1]+":01"
                            if len(type) == 3 and type[1] == "6" and type[2] == "8":
                                new_type = "A68:01"
                            if len(type) == 3 and type[1] == "2" and type[2] == "4":
                                new_type = "A24:01"
                            if len(type) == 3 and type[1] == "0" and type[2] == "2":
                                new_type = "A02:01"
                            if len(type) > 3 and "2" in new_type:
                                new_type = "A02:01"
                            if len(type) == 5 or (len(type) == 6 and type[5] == ","):
                                if type[2] == ":":
                                    new_type = type[0]+"0"+type[1]+":"+type[3]+type[4]
                                else:
                                    new_type = type[0]+type[1]+type[2]+":"+type[3]+type[4]
                            if len(type) == 6 and type[1] == "2":
                                new_type = "A02:01"

                            if ("(" in new_type or ")" in new_type) and "2" in new_type:
                                new_type = "A02:01"

                    new_key = new_key + new_type
                    pep = find_peptide_seq(filename)
                    ic50_value = 0
                    if pep != None:
                        ic50_value = get_ic50(new_key,pep)
                        result.append(filename+","+new_key+","+str(pep)+","+str(ic50_value))
                        #print(filename+","+org_key+","+new_key+","+str(pep)+","+str(ic50_value))
                    break
            filehandler.close()


    filehandler = open("ic50_HLA_A.csv", "w")
    for i in range(0, len(result)):
        filehandler.write(result[i]+"\n")
    filehandler.close()
    
