#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: Jan 16, 2017
#   Email: snerli@ucsc.edu
#

# bio libraries
from Bio.SubsMat.MatrixInfo import *

# import other required libraries
from collections import defaultdict


class BLOSUM:

    BLOSUMX = defaultdict(dict)

    def __init__(self, type=62):

        if type == 62:
            blosum = blosum62
        elif type == 80:
            blosum = blosum80
        elif type == 45:
            blosum = blosum45
        else:
            print("Unknown substitution matrix")

        for key in blosum:
            key0 = key[0]
            key1 = key[1]
            self.BLOSUMX[key0][key1] = blosum[key]
            self.BLOSUMX[key1][key0] = blosum[key]

        starArr = {}

        for key in self.BLOSUMX:
            starArr[key] = -4
            starArr['*'] = 1

        for key in starArr:
            self.BLOSUMX['*'][key] = starArr[key]
            self.BLOSUMX[key]['*'] = starArr[key]

    def get_matrix(self):
        return self.BLOSUMX
