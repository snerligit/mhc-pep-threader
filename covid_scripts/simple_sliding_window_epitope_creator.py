'''

Method to read multiple fasta files in a directory and create epitopes of given lengths

'''


import os
import sys
import argparse
from thread.fasta import FASTA

# read arguments
def get_args():
    """
    Parse command line arguments
    """

    parser = argparse.ArgumentParser(description="Method to read multiple fasta files and create epitope list of given lengths")
    parser.add_argument("-i", "--infile", help="fully qualified path to the fasta file")
    parser.add_argument("-l", "--length", type=int, default=9, help="Epitope length")

    args = parser.parse_args()

    return args

# perform sliding window and extract the epitopes of specified lengths
def sliding_window(args, fasta):

    headers = fasta.get_headers()
    epitopes = {}

    repeat = 0
    for header in headers:
        sequence = list(fasta.get_sequence(header))
        for i in range(0, len(sequence)-args.length+1):
            pep = ''.join(sequence[i:i+args.length])
            key = '>'+pep
            if key not in epitopes:
                epitopes[key] = pep
            else:
                repeat += 1

    for pep in epitopes:
        print (pep+"\n"+epitopes[pep])

    #print ("Total:", len(epitopes), "Repeats: ", repeat)

# main method
def main():

    args = get_args()

    fasta = FASTA(args.infile)
    fasta.read()

    sliding_window(args, fasta)


# Entry to the program
if __name__ == "__main__":
    main()
