#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: April 10, 2017
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

# import rosetta files
from thread.template import TEMPLATE
from thread.thread import THREAD
from thread.fasta import FASTA
from thread.pre_thread import PRE_THREADING
from input_output.input.mhc import MHC
from input_output.input.peptides import PEPTIDE
from input_output.input.beta2m import BETA2M
from input_output.input.chaperone import CHAPERONE
from input_output.input.tcr import TCR

# import other required libraries
import os
import sys
import imp

'''

MODEL class contains all the necessary functionalities required to
invoke necessary functionailities to carry out the p/MHC-TCR/Chaperone threading.

'''

class MODEL:

    # class members
    complex_headers = [] # headers of the complex
    complex_sequences = [] # sequences of the complex
    template = None # template structure
    args = None # input arguments
    mpi_install = True # flag enabled if mpi is installed

    # constructor
    def __init__(self, args):
        self.args = args

    # method to model mhcs for given peptide, beta2m, tcr or chaperone sequences
    def model_mhc_for_each_peptide_beta2m_tcr_chaperone(self):
        mhc = MHC(self.args.get_mhcs(), self.args.is_no_trim_mhc_flag_set())
        beta2m = BETA2M(self.args.get_beta2m())
        peptides = PEPTIDE(self.args.get_peptides())
        chaperones = CHAPERONE(self.args.get_chaperone())
        tcr = TCR(self.args.get_tcr())

        if ( len(peptides.get_headers()) > 1 ):
            if ( self.args.get_pep_start_index() == 0 ):
                print ("Please provide peptide start index. Option: -pep_start_index")
                exit(1)
            if ( self.args.get_interface_cupoint() == 0 ):
                print ("Please provide interface cutpoint. Option: -interface_cutpoint")
                exit(1)
        '''
            # Not sure why I introduced this. Keeping it commented until I figure out that case

            elif ( len(peptides.get_headers()) == 1 ):
            for pep_header in peptides.get_headers():
                if ( pep_header is not "none" ):
                    print ("Please provide peptide start index. Option: -pep_start_index")
                    exit(1)
        '''

        # get all combinations of sequences
        # mhcs, beta2ms, peptides, tcrs and chaperones
        for mhc_header in mhc.get_headers():
            for beta2m_header in beta2m.get_headers():
                for pep_header in peptides.get_headers():
                    for tcr_header in tcr.get_headers():
                        for chaperone_header in chaperones.get_headers():
                            header = mhc_header+"_"+beta2m_header+"_"+pep_header+"_"+tcr_header+"_"+chaperone_header
                            sequence = mhc.get_sequence(mhc_header)+beta2m.get_sequence(beta2m_header)+peptides.get_sequence(pep_header)+tcr.get_sequence(tcr_header)+chaperones.get_sequence(chaperone_header)
                            self.generate_fasta(header, sequence)

        # identify the number of jobs
        njobs = len(self.complex_sequences)

        try:
            from jd.job_distributor import JOB_DISTRIBUTOR
            job_dist = JOB_DISTRIBUTOR()
            # Can we make this call by only once processor?
            '''
            template_name = job_dist.perform_single_operation(self.get_template)
            if self.template == None and isinstance(template_name, str):
                self.template = TEMPLATE(template_name).get_pose_from_pdb()
                job_dist.apply(njobs, self.apply)
            elif self.template != None:
                job_dist.apply(njobs, self.apply)
            '''
            self.get_template()
            #job_dist.SimpleMPIJobDistributor(njobs, self.apply)
            job_dist.MPIJobDistributor(njobs, self.apply)
        except ImportError:
            self.mpi_install = False

        # if not mpi then run this on a single processor
        if not self.mpi_install:
            self.get_template()
            for job_id in range(njobs):
                self.apply(job_id)
                job_id += 1

    # apply method that calls threader
    def apply(self, i):
        pre_thread = None
        pre_thread = PRE_THREADING(self.template, self.complex_sequences[i],
                            self.complex_headers[i], self.args)
        threader = THREAD(pre_thread)
        threader.apply()

    # method to initialize template structure
    def get_template(self):
        self.template = TEMPLATE(self.args.get_template_pdb())
        self.template.treat_template_structure(self.args.get_mhc_chain(),
                                            self.args.get_peptide_chain(),
                                            self.args.is_no_trim_mhc_flag_set(),
                                            self.args.get_mhc_trim_length(),
                                            self.args.get_idealize_relax())
        return self.template.get_pdb()

    # method to generate fasta file names
    def generate_fasta(self, header, sequence):
        fasta_file_name = header+".fasta"
        FASTA(fasta_file_name).write(header, sequence)
        self.complex_headers.append(header)
        self.complex_sequences.append(sequence)
        return fasta_file_name
