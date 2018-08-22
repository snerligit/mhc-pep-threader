#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: August 22, 2018
#   Email: snerli@ucsc.edu
#

class HEADER_READER:

    filename = ""
    headers = []

    def __init__(self, filename):
        self.filename = filename
        self.headers = []

    def read_headers(self):
        if self.filename != None:
            read_file_handler = open(self.filename, "r")
            for line in read_file_handler:
                line = line.rstrip()
                self.headers.append(line)
            read_file_handler.close()

    def get_headers(self):
        if len(self.headers) == 0:
            self.headers.append('none')
        return self.headers
