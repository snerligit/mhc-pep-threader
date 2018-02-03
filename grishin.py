#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 26, 2017
#   Email: snerli@ucsc.edu
#

# additional bio libraries

# import other required libraries
import os
import sys
import subprocess

class GRISHIN:

    filename = ""
    target_head = ""
    template_head = ""
    target_string = ""
    template_string = ""

    def __init__(self, filename, target_head, template_head, target_string, template_string):
        self.filename = filename
        self.target_head = target_head
        self.template_head = template_head
        self.target_string = target_string
        self.template_string = template_string

    def get_file_name(self):
        return self.filename+".grishin"

    def write(self):
        writefile = open(self.get_file_name(), "w")
        writefile.write("## "+self.target_head+" "+self.template_head+"\n"+"#"+"\n")
        writefile.write("scores_from_program: 0\n")
        writefile.write("0 "+self.target_string+"\n")
        writefile.write("0 "+self.template_string+"\n")
        writefile.close()
