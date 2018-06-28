#/usr/bin/env python

import numpy as np
import logging
import time
import itertools
import sqlite3
import os
import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI

from utils.args import Args
from utils.logs import MyLogger
from utils import project
from utils import iotools
from qstats.jpa import JPA
from qstats.revreg import RevReg
from iotools.data import Data
from iotools.outhandler import Outhandler

# ==================================================
# Start MPI calculation
# =================================================


MPI.Init()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncore = comm.Get_size()

if rank == 0: start_time = time.time()

args = Args(comm, rank)
logger = MyLogger(__name__)

gtcent = None
gtnorm = None
expr   = None
maf = None

if rank == 0:
    data = Data(args)
    if args.simulate:
        data.simulate()
    else:
        data.load()
    gtcent = data.geno_centered
    logger.debug("Retained {:d} SNPs in {:d} samples".format(gtcent.shape[0], gtcent.shape[1]))
    gtnorm = data.geno_normed
    snpinfo = data.snpinfo
    expr = data.expression
    geneinfo = data.geneinfo
    maf = np.array([x.maf for x in snpinfo])
    
    #gtcent_shuf = np.zeros_like(gtcent)
    #gtnorm_shuf = np.zeros_like(gtnorm)
    #for i in range(gtcent.shape[0]):
    #    idx = np.random.permutation(np.arange(0,gtcent.shape[1]))
    #    np.random.shuffle(idx)
    #    gtcent_shuf[i,:] = gtcent[i,idx]
    #    gtnorm_shuf[i,:] = gtnorm[i,idx]    
    #gtcent = gtcent_shuf
    #gtnotm = gtnorm_shuf		################## SHUFFLED GENO ASSIGNMENT ALERT ! #######################

gtnorm = comm.bcast(gtnorm, root = 0)
gtcent = comm.bcast(gtcent, root = 0)
expr   = comm.bcast(expr,  root = 0)
maf  = comm.bcast(maf, root = 0)
comm.barrier()

if rank == 0: read_time = time.time()

jpa = JPA(gtnorm, expr, comm, rank, ncore, args.jpa)
jpa.compute()

if rank == 0: jpa_time = time.time()

# Select the SNPs with JPA score above threshold for RevReg
if args.jpa and args.rr:
    select = None
    if rank == 0:
        select = np.where(jpa.scores > args.jpacut)[0]
        logger.debug('JPA threshold: {:g}, and number of SNPs retained: {:d}'.format(args.jpacut, select.shape[0]))

    qselect = comm.bcast(select, root = 0)
    comm.barrier()
    #geno = geno[qselect, :]


if args.rr:
    sigbeta2 = np.repeat(args.sigmabeta ** 2, gtnorm.shape[0])
    if args.nullmodel == 'maf':
        rr = RevReg(gtnorm, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel, maf = maf)
    elif args.nullmodel == 'perm':
        #if rank == 0:
            #print(gtcent)
            #print(expr)
        rr = RevReg(gtcent, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel)
    rr.compute()

    
# Output handling only from master node // move it to module
if rank == 0: 
    if args.outprefix is None:
        args.outprefix = "out"    
    rr_time = time.time()
    if rank == 0:
        ohandle = Outhandler(args, snpinfo, geneinfo)
    if args.jpa:
        ohandle.write_jpa_out(jpa)
    if args.rr:
        ohandle.write_rr_out(jpa, rr)

if rank == 0: write_time = time.time()

if rank == 0:
    print ("File reading time: {:g} seconds".format(read_time - start_time))
    print ("JPA calculation time: {:g} seconds".format (jpa_time - read_time))
    print ("RR calculation time: {:g} seconds".format (rr_time - jpa_time))
    print ("Result writing time: {:g} seconds".format(write_time - rr_time))
    print ("Total execution time: {:g} seconds".format(time.time() - start_time))

MPI.Finalize()
