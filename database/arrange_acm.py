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
        head = "hla_sequences['"+line[1:]+"']="
    else:
        index = line.find('SHSM')
        if index >= 0 and len(line) >300:
            head = head+"'G"+line[index:index+275]+"'"
        else:
            pass

        if head.find('GSHS') != -1:
            print head

fastafilehandler.close()
