#/usr/bin/env python

import numpy as np
import logging
import time
import itertools
import os
import gc
import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI

from utils.args import Args
from utils.logs import MyLogger
from utils import project
from utils import mpihelper
from iotools.data import Data

from qstats.zstats import ZSTATS
from qstats.nullscores import JPANULL
from scipy.linalg import eigh

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
snpinfo = None
masklist = None
maskcomp = None
W = None
Q = None
Zmean = None
niter = None

if rank == 0:
    logger.debug("Using {:d} cores".format(ncore))
    data = Data(args)
    if args.simulate:
        data.simulate()
    else:
        data.load()
    tissue = args.gx_file.split(".")[4]
    logger.debug("Tissue: "+tissue)
    gtcent = data.geno_centered
    gtnorm = data.geno_normed
    outdir = os.path.dirname(args.outprefix)
    maxmemSNPs = 20000000 #100000000
    # maxmemSNPs = 1440000 #360*4000
    nchunks = 1
    if ((gtcent.shape[0] * gtcent.shape[1]) > maxmemSNPs): # is a conservative number assuming np.float64 uses 8bytes ~3GB object
        # calculate how many chunks of SNPs we will need
        print(gtnorm.shape)
        chunk_size = int(maxmemSNPs / gtnorm.shape[1])
        nchunks = int(gtnorm.shape[0] / chunk_size ) + 1
        print("Going to iterate over {:g} chunks of size {:g} SNPs".format(nchunks, chunk_size))
        gt_chunks, offset = mpihelper.split_genotype(gtnorm, nchunks)
        print(len(gt_chunks))
        print(offset)
        print([x.shape for x in gt_chunks])
    snpinfo = data.snpinfo
    expr = data.expression
    geneinfo = data.geneinfo
    logger.debug("After prefilter: {:d} SNPs and {:d} genes in {:d} samples".format(gtcent.shape[0], expr.shape[0], gtcent.shape[1]))
#<<<<<<< Updated upstream
    niter = 100000
else:
    gtnorm = None
    expr   = None
    gt_chunks = None
    nchunks = None

expr    = comm.bcast(expr,  root = 0)
nchunks = comm.bcast(nchunks, root = 0)
if rank != 0:
    gt_chunks = [None for i in range(nchunks)]
#=======
#    niter = 10000
#
#gtnorm = comm.bcast(gtnorm, root = 0)
#gtcent = comm.bcast(gtcent, root = 0)
#expr   = comm.bcast(expr,  root = 0)
#>>>>>>> Stashed changes

snpinfo = comm.bcast(snpinfo, root = 0)
comm.barrier()

if rank == 0: read_time = time.time()
if rank == 0: logger.debug("Computing Z-stats")

if nchunks > 1:
    zstats_buf = None
    for i in range(nchunks):
        if rank == 0: print(len(gt_chunks[i]))
        zstats = ZSTATS(gt_chunks[i], expr, comm, rank, ncore)
        zstats.compute()

        if rank == 0:
            np.savetxt(outdir+"/chunk{0:03d}_Zstats.txt".format(i), zstats.scores)
            if zstats_buf is None:
                zstats_buf = zstats.scores
            else:
                zstats_buf = np.vstack((zstats_buf, zstats.scores))
            logger.debug("Zstats: "+str(zstats_buf.shape))
else:
    zstats = ZSTATS(gtnorm, expr, comm, rank, ncore, masklist)
    zstats.compute()
    if rank == 0: zstats_buf = zstats.scores

if rank == 0: zstat_time = time.time()
if rank == 0: logger.debug("Computing W and Q")

if rank == 0:
    logger.debug("Final Zstats: "+str(zstats_buf.shape))
    C = np.cov(zstats_buf.T)
    # Numpy gives imaginary eigenvalues, use eigh from scipy
    # for decomposition of real symmetric matrix
    W, Q = eigh(C)
    logger.info("W and Q calculation took: {:g} seconds".format(time.time() - zstat_time))
    logger.info("W: {:s}".format(str(W.shape)))
    logger.info("Q: {:s}".format(str(Q.shape)))
    # still some eigenvalues are negative. force them to zero if they are negligible. (!!!!!!!!!!!)
    # check if everything is ok
    Wsparse = W.copy()
    Wsparse[np.where(W < 0)] = 0
    
    if not np.allclose(C, Q @ np.diag(Wsparse) @ Q.T):
        logger.error("Eigen vectors could not be forced to positive")
        exit
    else:
        W = Wsparse
        logger.debug("Eigen vectors are forced to positive")
        np.savetxt(outdir+"/Wcpma_"+tissue+".txt", W)
        np.savetxt(outdir+"/Qcpma_"+tissue+".txt", Q)
    Zmean = np.mean(zstats, axis = 0)

comm.barrier()

W = comm.bcast(W, root = 0)
Q = comm.bcast(Q, root = 0)
Zmean = comm.bcast(Zmean, root = 0)
niter = comm.bcast(niter, root = 0)
comm.barrier()

raise

if rank == 0: wq_time = time.time()
if rank == 0: logger.debug("Computing null JPA-scores")

jpa_null = JPANULL(W, Q, Zmean, niter, comm, rank, ncore)
jpa_null.compute()

if rank == 0: null_time = time.time()
if rank == 0:
    fname = args.qnullfile
    with open(fname, 'w') as fout:
        for qnull in jpa_null.scores:
            fout.write("{:g}\n".format(qnull)) 

if rank == 0: write_time = time.time()

if rank == 0:
    logger.info("File reading time: {:g} seconds".format(read_time - start_time))
    logger.info("F-stat calculation time: {:g} seconds".format (zstat_time - read_time))
    logger.info("W and Q calculation time: {:g} seconds".format(wq_time - zstat_time))
    logger.info("Null JPA calculation time: {:g} seconds".format (null_time - wq_time))
    logger.info("Result writing time: {:g} seconds".format(write_time - null_time))
    logger.info("Total execution time: {:g} seconds".format(time.time() - start_time))

MPI.Finalize()
