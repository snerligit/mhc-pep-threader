import os
import sys
import numpy
import argparse
from collections import defaultdict

'''
python ~/rosettamhc/mhc-pep-threader/covid/identify_common_cold_peps.py -e sars_eptiopes.txt

'''

# get command line arguments
def get_args():
    """
    Parse command line arguments
    """

    parser = argparse.ArgumentParser(description="Method to setup rosettamhc jobs")
    parser.add_argument("-e", help="sars eptiopes, one line for each epitope in the text file")
    parser.add_argument("-aln", help="alignment file")
    parser.add_argument("-tag", help="protein name")
    parser.add_argument("-length", help="epitope length", default=9, type=int)
    parser.add_argument("-name", help="protein name [229E, NL63, OC43, HKU1]")

    args = parser.parse_args()
    return args

# read SARS-CoV-2 epitope file
def read_epitopes(args):

    inputfilehandler = open(args.e, 'r')
    epitopes = []
    for line in inputfilehandler:
        line = line.rstrip()
        epitopes.append(line)
    inputfilehandler.close()

    return epitopes

# perform sliding window
def sliding_window(seq, args):

    sequence = list(seq)
    epitopes = {}
    repeat = 0
    length = args.length

    for i in range(0, len(sequence)-length+1):
        pep = ''.join(sequence[i:i+length])
        key = pep
        if key not in epitopes:
            epitopes[pep] = i
        else:
            repeat += 1

        #print ("Found repeats: ", repeat)

    return epitopes

# read alignment file
def read_aln(args):

    inputfilehandler = open(args.aln, 'r')
    sarscov2 = ''
    hcov229e = ''
    hcovnl63 = ''
    hcovoc43 = ''
    hcovhku1 = ''
    for line in inputfilehandler:
        line = line.rstrip()
        if "CLUSTAL" not in line and line != '' and '*' not in line and ':' not in line:
            fields = line.split()
            if '229E_strain' in line:
                hcov229e += fields[1]
            elif 'NL63_strain' in line:
                hcovnl63 += fields[1]
            elif 'SARS_COV2' in line:
                sarscov2 += fields[1]
            elif 'OC43_strain' in line:
                hcovoc43 += fields[1]
            elif 'HKU1_strain' in line:
                hcovhku1 += fields[1]
    inputfilehandler.close()

    return (sarscov2, hcov229e, hcovnl63, hcovoc43, hcovhku1)

def print_epitopes(e, match, args, strain):

    if '-' not in match:
        print ('>'+e+'_'+match+'_'+args.tag+'_'+strain)
        print (match)

def match(epitopes, new_epitopes, hcov229e, hcovnl63, hcovoc43, hcovhku1, args):

    hcov229e_arr = list(hcov229e)
    hcovnl63_arr = list(hcovnl63)
    hcovoc43_arr = list(hcovoc43)
    hcovhku1_arr = list(hcovhku1)

    specific = {}

    print ("SARS-CoV-2 cross-reactive epitopes:")
    for e in epitopes:
        e_len = len(list(e))
        if e in new_epitopes:
            hcov229e_match = ''.join(hcov229e_arr[new_epitopes[e]:new_epitopes[e]+e_len])
            hcovnl63_match = ''.join(hcovnl63_arr[new_epitopes[e]:new_epitopes[e]+e_len])
            hcovoc43_match = ''.join(hcovoc43_arr[new_epitopes[e]:new_epitopes[e]+e_len])
            hcovhku1_match = ''.join(hcovhku1_arr[new_epitopes[e]:new_epitopes[e]+e_len])

            if args.name == '229E':
                print_epitopes(e, hcov229e_match, args, "229E")
            elif args.name == 'NL63':
                print_epitopes(e, hcovnl63_match, args, "NL63")
            elif args.name == 'OC43':
                print_epitopes(e, hcovoc43_match, args, "OC43")
            elif args.name == 'HKU1':
                print_epitopes(e, hcovhku1_match, args, "HKU1")

            if '-' in hcov229e_match or '-' in hcovhku1_match or '-' in hcovnl63_match or '-' in hcovoc43_match:
                specific[e] = e+','+args.tag+','+hcov229e_match+','+hcovhku1_match+','+hcovnl63_match+','+hcovoc43_match
            #print ("Found "+e, hcov229e_match, hcovnl63_match, hcovoc43_match, hcovhku1_match)

    print ("SARS-CoV-2 specific epitopes:")
    for pep in specific:
        print (specific[pep])

# main method
def main(args):

    epitopes = read_epitopes(args)
    (sarscov2, hcov229e, hcovnl63, hcovoc43, hcovhku1) = read_aln(args)
    new_epitopes = sliding_window(sarscov2, args)
    match(epitopes, new_epitopes, hcov229e, hcovnl63, hcovoc43, hcovhku1, args)

# Entry to the program
if __name__ == "__main__":

    main(get_args())
