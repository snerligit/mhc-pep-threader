#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: December 1, 2017
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

# import rosetta files
from pyrosetta.rosetta.protocols.comparative_modeling import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.core.sequence import *
from pyrosetta.rosetta.protocols.jobdist import *

# additional bio libraries

#custom libraries
from relax import RELAX
from template import TEMPLATE
from fasta import FASTA
from grishin import GRISHIN
from align import ALIGN
from movemap import MOVEMAP
from job_distributor import JOB_DISTRIBUTOR
from chain_split import CHAIN_SPLIT
from interface_analyzer import INTERFACE
from job_distributor import JOB_DISTRIBUTOR

# import other required libraries
import os
import sys


class THREAD:

    target = None
    template = None
    alignment = None
    grishin_file_name = None
    fasta = ""
    nstruct = 1
    peptide = ""
    pep_header = ""
    mhc_headers = []
    #pep_start_index = 181
    #pep_start_index = 376
    pep_start_index = 181
    groove_distance = 3.5
    pep_length = 0

    def __init__(self, template, fasta, nstruct, pep_header, peptide, pep_start_index):

        self.template = template
        self.fasta = fasta
        self.nstruct = nstruct
        self.peptide = peptide
        self.pep_header = pep_header
        #self.pep_start_index = 181
        self.pep_start_index = pep_start_index
        self.groove_distance = 3.5
        self.pep_length = len(list(peptide))

    def apply(self):

        self.align_template_target_sequences()
        self.do_homology_modeling()

    def get_target_file_name(self, header="HLAs"):

            return header+"_on_"+self.template.get_stripped_name()+"_with_"+self.pep_header

    def get_pep_start_index(self, sequence):

        index = 0
        print (sequence)
        print (self.peptide)
        try:
            index = sequence.rindex(self.peptide)
            return index
        except ValueError:
            print("Retry with different alignment scheme")
            exit(1)

    def align_template_target_sequences(self, clustal=True):

        self.target = FASTA(self.fasta)

        if clustal == False:
            self.target.read()
            target_seqs = self.target.get_sequences()

            for key,value in target_seqs.items():

                self.alignment = ALIGN(self.template.get_sequence(), self.target.get_sequence(key))
                self.alignment.align()

                # template_seq = query_sequence , target_seq =  target_sequence
                grishin = GRISHIN(self.get_target_file_name(), key, self.template.get_name(),
                        self.alignment.get_aligned_target(), self.alignment.get_aligned_query())

                self.mhc_headers.append(key)
                self.grishin_file_name = grishin.get_file_name()
                grishin.write()
        else:
            self.alignment = ALIGN(self.template.get_sequence(), None, self.fasta)
            self.alignment.clustal()
            alignment = FASTA(self.alignment.get_clustal_output_filename())
            alignment.read()
            aln_seqs = alignment.get_sequences()

            for key,value in aln_seqs.items():
                if key != 'template':
                    grishin = GRISHIN(self.get_target_file_name(), key, self.template.get_name(),
                             value, aln_seqs['template'])

                    self.mhc_headers.append(key)
                    self.grishin_file_name = grishin.get_file_name()
                    grishin.write()

    def do_homology_modeling(self):

        #jd = JOB_DISTRIBUTOR(self.get_target_file_name())
        #jd.set_native(self.template.get_pose())
        #silent_file = SILENT()
        alignment_vec = read_aln("grishin", self.grishin_file_name)
        ctr = 0
        for vec in alignment_vec:
            print(vec)
            new_pose = pose_from_sequence(self.target.get_sequence(self.mhc_headers[ctr]), "fa_standard")
            thread = PartialThreadingMover(vec, self.template.get_pose())
            thread.apply(new_pose)

            #chain = CHAIN_SPLIT(new_pose, self.get_pep_start_index(new_pose.sequence()))
            #print (self.get_pep_start_index(new_pose.sequence()))
            #chain.cut()
            #threaded = chain.get_pose()
            threaded = new_pose
            tag = self.get_target_file_name(self.mhc_headers[ctr])
            #jd.get_dist_obj.output_decoy(threaded)
            threaded.dump_pdb(tag+".pdb")
            movemap = MOVEMAP(tag+".pdb", self.pep_start_index, self.pep_length, self.groove_distance, tag+".movemap", True)
            movemap.apply()
            #silent_file.add(threaded,tag,"hla_output.out")
            ctr += 1
