
# Author: Santrupti Nerli
# Date: Feb 2017

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

align_file = "aligned_groove_3.5.clustal"
identity_file = "sequence_identity_B15107.csv"

readFile = open(align_file, "r")
writeFile = open(identity_file, "w")

target_seq = {}
target_name = {}

i = 0
temp = -1
template_seq = ""
for line in readFile:
    if i == 3105:
        i = 0
    modified = re.sub(r"\s+", " ", line)
    pair = modified.split(" ")
    #if("B_15_107" in line): # 9-mer

    # this is where you need to provide the pattern of template file name
    # For example, I had B_15_84 pattern for the template file. You can see this in the
    # aligned_groove_3.5.clustal file
    if ("B_15_84" in line): # 10-mer
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
        name = re.sub(".pdb", "", target_name[t])

        # The below lines are used because I used a template_A_01_01 naming convention
        # and I had to convert them to A*01:01 convention
        # you can adjust the following names to whatever convention you want

        if ("template_A" in name) or ("template_B" in name) or ("template_C" in name):
            name = re.sub(r"template_", "", name)
            nameArr = name.split("_")
            name = nameArr[0]+"*"+nameArr[1]+":"+nameArr[2]
            #print name
        else:
            name = re.sub(r"template_", "B*15:", name)
            arr = name.split("_")
            name = arr[0]+"*"+arr[1]+":"+arr[2]
        if(temp == t):
            print(score)
        writeFile.write(name+","+str(score)+"\n")

writeFile.close()
readFile.close()
