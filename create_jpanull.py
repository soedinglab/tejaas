#/usr/bin/env python

import numpy as np
import logging
import time
import itertools
import os
import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI

from utils.args import Args
from utils.logs import MyLogger
from utils import project
from iotools.data import Data
#from iotools.outhandler import Outhandler
#from iotools import readmaf

from qstats.jpa_zstats import JPAZSTATS
from qstats.jpa_null import JPANULL

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
#maf = None
snpinfo = None
#masklist = None
#maskcomp = None
W = None
Q = None
niter = None

if rank == 0:
    logger.debug("Using {:d} cores".format(ncore))
    data = Data(args)
    if args.simulate:
        data.simulate()
    else:
        data.load()
    gtcent = data.geno_centered
    gtnorm = data.geno_normed
    snpinfo = data.snpinfo
    expr = data.expression
    geneinfo = data.geneinfo
    #masklist = data.cismasks_list
    #maskcomp = data.cismasks_comp
    logger.debug("After prefilter: {:d} SNPs and {:d} genes in {:d} samples".format(gtcent.shape[0], expr.shape[0], gtcent.shape[1]))
    
    #maf = readmaf.load(snpinfo, args.nullmodel, args.maf_file)
    niter = 1000

gtnorm = comm.bcast(gtnorm, root = 0)
gtcent = comm.bcast(gtcent, root = 0)
expr   = comm.bcast(expr,  root = 0)
snpinfo = comm.bcast(snpinfo, root = 0)
#maf  = comm.bcast(maf, root = 0)
#masklist = comm.bcast(masklist, root = 0)
#maskcomp = comm.bcast(maskcomp, root = 0)
comm.barrier()

if rank == 0: read_time = time.time()
if rank == 0: logger.debug("Computing Z-stats")

jpa_zstats = JPAZSTATS(gtnorm, expr, comm, rank, ncore)
jpa_zstats.compute()

if rank == 0: zstat_time = time.time()
if rank == 0: logger.debug("Computing W and Q")

if rank == 0:
    zstats = jpa_zstats.zstats
    C = np.cov(zstats.T)
    W, Q = np.linalg.eig(C)

W = comm.bcast(W, root = 0)
Q = comm.bcast(Q, root = 0)
niter = comm.bcast(niter, root = 0)
comm.barrier()

if rank == 0: wq_time = time.time()
if rank == 0: logger.debug("Computing null JPA-scores")

jpa_null = JPANULL(W, Q, niter, comm, rank, ncore)
jpa_null.compute()

if rank == 0: null_time = time.time()
if rank == 0:
    fname = args.qnullfile
    with open(fname, 'w') as fout:
        for qnull in jpa_null.scores:
            fout.write(f"{qnull}\n") 

if rank == 0: write_time = time.time()

if rank == 0:
    logger.info("File reading time: {:g} seconds".format(read_time - start_time))
    logger.info("F-stat calculation time: {:g} seconds".format (zstat_time - read_time))
    logger.info("W and Q calculation time: {:g} seconds".format(wq_time - zstat_time))
    logger.info("Null JPA calculation time: {:g} seconds".format (null_time - wq_time))
    logger.info("Result writing time: {:g} seconds".format(write_time - null_time))
    logger.info("Total execution time: {:g} seconds".format(time.time() - start_time))

MPI.Finalize()
