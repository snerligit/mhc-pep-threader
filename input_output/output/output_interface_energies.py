#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: February 2, 2018
#   Email: snerli@ucsc.edu
#

'''

OUT_INTERFACE_ENERGY class contains all the necessary functionalities required
create output file containing targets and their binding energies.

'''

# import required libraries
from collections import defaultdict

class OUT_INTERFACE_ENERGY:

    # class members
    filename = "" # filename to write the output to
    energies = None # dictionary to store the binding energies

    # constructor
    def __init__(self, filename):
        self.energies = defaultdict(list)
        self.filename =  filename

    # method to populate dictionary with binding energies
    def add(self, key, value, rms1 = -1, rms2 = -1):
        fields = key.split("_")[:4]
        if rms1 == -1 and rms2 == -1:
            self.energies[key].append(value)
        else:
            metrics = []
            metrics.append(value)
            metrics.append(rms1)
            metrics.append(rms2)
            self.energies[key].append(metrics)

    # method to output a csv file containing names of domains
    # and their corresponding binding energies
    def write(self):
        write_file_handler = open(self.filename, "a")
        for key,value in self.energies.items():
            fields = key.split("_")
            mhc = fields[0]
            beta2m = fields[1]
            pep = fields[2]
            tcr = fields[3]
            chaperone = fields[4]
            energy = ""
            for i in range(len(value)):
                energy += str(value[i])
                if i != len(value)-1:
                    energy += ","
            text = mhc+","+beta2m+","+pep+","+tcr+","+chaperone+","+energy+"\n"
            write_file_handler.write(text)
            print (text)
        write_file_handler.close()
