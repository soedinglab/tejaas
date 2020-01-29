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
from qstats.jpa import JPA
from qstats.revreg import RevReg
from qstats.jpa_pvals import JPAPVALS

from iotools.data import Data
from iotools.data import optimize_sb2
from iotools.outhandler import Outhandler
from iotools import readmaf
from iotools import readqnull

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
snpinfo = None
masklist = None
maskcomp = None
if args.onlyjpa: qnull = None

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
    masklist = data.cismasks_list
    maskcomp = data.cismasks_comp
    logger.debug("After prefilter: {:d} SNPs and {:d} genes in {:d} samples".format(gtcent.shape[0], expr.shape[0], gtcent.shape[1]))
    
    maf = readmaf.load(snpinfo, args.nullmodel, args.maf_file)
    if args.onlyjpa: qnull = readqnull.load(args.qnullfile)

gtnorm = comm.bcast(gtnorm, root = 0)
gtcent = comm.bcast(gtcent, root = 0)
expr   = comm.bcast(expr,  root = 0)
snpinfo = comm.bcast(snpinfo, root = 0)
maf  = comm.bcast(maf, root = 0)
masklist = comm.bcast(masklist, root = 0)
maskcomp = comm.bcast(maskcomp, root = 0)
if args.onlyjpa: qnull = comm.bcast(qnull, root = 0)
comm.barrier()

if rank == 0: read_time = time.time()

if rank == 0: logger.debug("Computing JPA")

if args.onlyjpa:
    jpa = JPAPVALS(gtnorm, expr, qnull, comm, rank, ncore, args.onlyjpa, masklist)
    jpa.compute()
else:
    jpa = JPA(gtnorm, expr, comm, rank, ncore, args.jpa, masklist)
    jpa.compute()

if rank == 0: jpa_time = time.time()

# # Select the SNPs with JPA score above threshold for RevReg
if args.jpa and args.rr:
    select = None
    if rank == 0:
        select = np.where(jpa.scores > args.jpacut)[0]
        logger.debug('JPA threshold: {:g}, and number of SNPs retained: {:d}'.format(args.jpacut, select.shape[0]))
    qselect = comm.bcast(select, root = 0)
    comm.barrier()
    #geno = geno[qselect, :]
    #broadcast this new genotype

if args.rr:
    if args.dynamic:
        if rank == 0:
            logger.debug("-- DINAMICALLY ADJUSTING SIGMA BETA | Target: {:f}-- ".format(args.dynamic))
        # adjust sigma_beta2 for each SNP (aprox)
        Yt = expr.T
        U, S, Vt  = np.linalg.svd(Yt, full_matrices=False)
        sigmax2   = np.var(gtcent, axis = 1)
        sigbeta2 = optimize_sb2(S, sigmax2, args.dynamic)

        _S2 = np.square(S)
        _S2mod = _S2 + (sigmax2[0] / sigbeta2[0])
        Keff = np.sum(_S2/_S2mod) / len(_S2)
        if rank == 0:
            logger.debug("Current Keff @ {:f}".format(Keff))
    elif args.mml:
        if rank == 0:
            logger.debug("SIGMA BETA OPTIMIZED BY MML")
        sigbeta2 = [None for x in range(gtnorm.shape[0])]
    else:
        if rank == 0:
            logger.debug("SIGMA BETA FIXED TO {:g}".format(args.sigmabeta))
        sigbeta2 = np.repeat(args.sigmabeta ** 2, gtnorm.shape[0])

    if args.nullmodel == 'maf':
        rr = RevReg(gtnorm, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel, maf = maf, masks = maskcomp)
    elif args.nullmodel == 'perm':
        rr = RevReg(gtcent, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel, maf = maf, masks = maskcomp)

    rr.compute(get_betas = False)

if rank == 0: rr_time = time.time()

# Output handling only from master node // move it to module
if rank == 0: 
    ohandle = Outhandler(args, snpinfo, geneinfo)
    if args.onlyjpa:
        ohandle.write_jpa_pvals(jpa)
    if args.rr:
        ohandle.write_rr_out(jpa, rr, write_betas = True)

if rank == 0: write_time = time.time()

if rank == 0:
    logger.info("File reading time: {:g} seconds".format(read_time - start_time))
    logger.info("JPA calculation time: {:g} seconds".format (jpa_time - read_time))
    logger.info("RR calculation time: {:g} seconds".format (rr_time - jpa_time))
    logger.info("Result writing time: {:g} seconds".format(write_time - rr_time))
    logger.info("Total execution time: {:g} seconds".format(time.time() - start_time))

MPI.Finalize()
