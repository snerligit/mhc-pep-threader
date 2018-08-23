#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 26, 2017
#   Email: snerli@ucsc.edu
#

'''

GRISHIN class contains all the necessary functionalities required to create
Rosetta specific alignment file.

'''

# import other required libraries
import os
import sys
import subprocess

class GRISHIN:

    # class members
    filename = "" # grishin file name
    target_head = "" # target name
    template_head = "" # template name
    target_sequence = "" # target sequence
    template_sequence = "" # template sequence

    # constructor
    def __init__(self, filename, target_head, template_head, target_sequence, template_sequence):
        self.filename = filename
        self.target_head = target_head
        self.template_head = template_head
        self.target_sequence = target_sequence
        self.template_sequence = template_sequence

    # method to create and write to the grishin file
    # the formatting is very specific to Rosetta
    # See the link: https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/Grishan-format-alignment
    def write(self, is_new = False):
        if is_new:
            writefile = open(self.get_file_name(), "w")
        else:
            writefile = open(self.get_file_name(), "a")
        writefile.write("## "+self.target_head+" "+self.template_head+"\n"+"#"+"\n")
        writefile.write("scores_from_program: 0\n")
        writefile.write("0 "+self.target_sequence+"\n")
        writefile.write("0 "+self.template_sequence+"\n")
        writefile.write("--\n")
        writefile.close()

    # getter method
    def get_file_name(self):
        return self.filename+".grishin"
