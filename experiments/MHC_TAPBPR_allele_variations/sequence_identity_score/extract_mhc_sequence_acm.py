import os
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description='arrange sequences')
parser.add_argument("-fasta", help="fasta file", required='True')

args = parser.parse_args()

fasta_ = args.fasta

fastafilehandler = open(fasta_, "r")

for line in fastafilehandler:
    line = line.rstrip()
    if ">" in line:
        head = ">mhcsequence_"+line[1:]+"\n"
    else:
        index = line.find('')
        if index >= 0 and len(line) >300:
            head = head+""+line[index:index+275]+""
        else:
            pass

        if head.find('') != -1:
            print head

fastafilehandler.close()
