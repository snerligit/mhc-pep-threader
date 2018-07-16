
# Author: Santrupti Nerli
# Date: Feb 2017
# Edit by ACM, June 2018

import re
import Bio
from Bio.SubsMat import MatrixInfo

def score_match(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]

def score_pairwise(seq1, seq2, matrix, gap_s, gap_e):
    score = 0
    gap = False
    for i in range(len(seq1)):
        print (seq2, i)
        pair = (seq1[i], seq2[i])
        if not gap:
            if '-' in pair:
                gap = True
                score += gap_s
            else:
                score += score_match(pair, matrix)
        else:
            if '-' not in pair:
                gap = False
                score += score_match(pair, matrix)
            else:
                score += gap_e
    return score

# Provide clustalomega alignment file.
align_file = "clustalo-I20180624-001744-0615-40632901-p1m.clustal"
# Provide name for the output file.
identity_file = "sequence_identity_A3001.csv"

readFile = open(align_file, "r")
writeFile = open(identity_file, "w")

target_seq = {}
target_name = {}

i = 0
temp = -1
template_seq = ""
for line in readFile:
    if i == 84:
        i = 0
    modified = re.sub(r"\s+", " ", line)
    pair = modified.split(" ")

    # Provide pattern for the template to align to.
    if ("A*30:01" in line):
        temp = i
        template_seq = template_seq+pair[1]
    if i not in target_name:
        target_name[i] = pair[0]
    if i in target_seq:
        target_seq[i] = target_seq[i]+pair[1]
    else:
        target_seq[i] = pair[1]
    i += 1

blosum = MatrixInfo.blosum62
for t in target_seq:
    if target_seq[t].rstrip() != "":
        score = score_pairwise(template_seq, target_seq[t], blosum, 0, 0)
        #name = re.sub(target_name[t])
        name = target_name[t]
        if(temp == t):
            print(score)
        writeFile.write(name+","+str(score)+"\n")

writeFile.close()
readFile.close()
