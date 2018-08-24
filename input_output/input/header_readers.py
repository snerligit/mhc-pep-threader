#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: August 22, 2018
#   Email: snerli@ucsc.edu
#

'''

HEADER_READER class contains all the necessary functionalities required to
read respective headers from the user input files.

'''

class HEADER_READER:

    # class members
    filename = "" # input file name containing beta2m, chaperone or mhc header lists
    headers = [] # header list

    # constructor
    def __init__(self, filename):
        self.filename = filename
        self.headers = []

    # method to read the file and populate headers
    # if filenames are not specified. For instance, chaperone,
    # then do nothing
    def read_headers(self):
        if self.filename != None:
            read_file_handler = open(self.filename, "r")
            for line in read_file_handler:
                line = line.rstrip()
                self.headers.append(line)
            read_file_handler.close()

    # getter method
    # if filename is not provided, then just return none
    def get_headers(self):
        if len(self.headers) == 0:
            self.headers.append('none')
        return self.headers
