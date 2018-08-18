#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 31, 2018
#   Email: snerli@ucsc.edu
#

# Load mpi4py library
from mpi4py import MPI

class JOB_DISTRIBUTOR:

    rank = -1

    def apply(self, njobs, fun):

        comm = MPI.COMM_WORLD

        print ("No. of jobs received: ", njobs )

        self.rank = comm.Get_rank()
        size = comm.Get_size()

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
            if i%size != self.rank: continue
            print ("Performing task ", i, " now in rank: ", self.rank )
            fun(i)
