#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 31, 2018
#   Email: snerli@ucsc.edu
#

# Load mpi4py library
from mpi4py import MPI
import pickle

class JOB_DISTRIBUTOR:

    rank = -1
    comm = None
    size = -1
    done = False

    def __init__(self):

        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

    def get_communicator(self):

        return self.comm

    def get_rank(self):

        return self.rank

    def sync(self):

        self.comm.barrier()

    def broadcast(self, data):

        return self.comm.bcast(data, root = 0)

    def set_all_dones(self):

        self.done = self.comm.bcast(self.done, root = 0)

    def perform_single_operation(self, fun):

        if not self.done:
            self.sync()
            ret_value = ""
            if self.rank == 0:
                print ("In rank: ", self.rank, " .Performing template preparation.")
                ret_value = fun()
                self.done = True
                self.set_all_dones()

            # broadcast the ret_val to other ranks
            ret_value = self.comm.bcast(ret_value, root=0)
            print ("Reached here. Return value is: ", ret_value)
            return ret_value


    def apply(self, njobs, fun):

        print ("No. of jobs received: ", njobs )

        '''
        myjobs = []

        if rank == 0:
            jobs = list(range(njobs))
            jobs.extend( [None]*(size - njobs % size) )
            n = len(jobs)/size
            print ("Size of n: ", n, " and jobs: ", len(jobs) )
            for i in range(size):
                queue = []  # list of jobs for individual cpu
                for j in range(int(n)):
                    print ("math: ", j*size+i)
                    queue.append(jobs[j*size+i])

                if( i == 0 ):
                    myjobs = queue
                else:
                    # now sending the queue to the process
                    logger.info('Sending %s to node %s' % (queue, i) )
                    comm.send(queue, dest=i)
        else:
            # getting decoy lists
            myjobs = comm.recv(source=0)

        print('Node %s, got queue:%s' % (rank, myjobs) )

        for j in myjobs:
            if j is not None and j < njobs:
                print (" the value of j is: ", j)
                fun(j)
        '''
        for i in range(njobs):
            if i%self.size != self.rank: continue
            print ("Performing task ", i, " now in rank: ", self.rank )
            fun(i)
