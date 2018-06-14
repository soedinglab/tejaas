#/usr/bin/env python

import numpy as np
import logging
import time

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

# ==================================================
# Start MPI calculation
# =================================================

start_time = time.time()

MPI.Init()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncore = comm.Get_size()

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
    gtnorm = data.geno_normed
    snpinfo = data.snpinfo
    expr = data.expression
    geneinfo = data.geneinfo
    maf = np.array([x.maf for x in snpinfo])

gtnorm = comm.bcast(gtnorm, root = 0)
gtcent = comm.bcast(gtcent, root = 0)
expr   = comm.bcast(expr,  root = 0)
maf  = comm.bcast(maf, root = 0)
comm.barrier()

jpa = JPA(gtnorm, expr, comm, rank, ncore, args.jpa)
jpa.compute()

# Select the SNPs with JPA score above threshold for RevReg
if args.jpa and args.rr:
    select = None
    if rank == 0:
        select = np.where(jpa.scores > args.jpacut)[0]
        logger.debug('JPA threshold: {:g}, and number of SNPs retained: {:d}'.format(args.jpacut, select.shape[0]))

    qselect = comm.bcast(select, root = 0)
    comm.barrier()
    geno = geno[qselect, :]

if args.rr:
    sigbeta2 = np.repeat(args.sigmabeta ** 2, gtnorm.shape[0])
    if args.nullmodel == 'maf':
        rr = RevReg(gtnorm, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel, maf = maf)
    elif args.nullmodel == 'perm':
        rr = RevReg(gtcent, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel)
    rr.compute()

# Output handling only from master node // move it to module
if rank == 0:
    pvals = jpa.pvals
    logger.debug('Pval matrix size: {:d} x {:d}'.format(pvals.shape[0], pvals.shape[1]))
#    outhandler = OutputHandler(outprefix)
#    outhandler.writedb(snpinfo, geneinfo, pvals)
    if args.jpa:
        jpascores = jpa.scores
#        outhandler.writejpa(snpinfo, jpascores)
#        logger.debug('Scores size: {:d}'.format(jpascores.shape[0]))
    if args.rr:
        rrscores = rr.scores
        pvals = rr.pvals
        mu = np.mean(rr.null_mu)
        sigma = np.mean(rr.null_sigma)
        logger.debug('Mean of RR scores: {:g}, Mean of RR null: {:g}\n'.format(np.mean(rrscores), mu))
        logger.debug('Variance of RR scores: {:g}, Variance of RR null: {:g}\n'.format(np.std(rrscores), sigma))
#        outhandler.writerr(snpinfo, rrscores)

if rank == 0:
    print ("Job completed in {:g} seconds".format(time.time() - start_time))

MPI.Finalize()
