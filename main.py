#/usr/bin/env python
'''
Tejaas: Discover Trans-eQTLs!

This initiates the MPI calls for Tejaas.
It is a wrapper for all the options used in Tejaas.
Currently, the following arguments are processed:
1) args.jpa
    This is very similar to the CPMA-score of Brynedal et. al.
    Requires a file containing null JPA-scores.
    The user must provide a filepath using the command: --jpanull filepath
        - If file exists, then it uses the file.
        - If file does not exist, then it creates a new null file at the specified filepath.
    Performs JPA-score analysis and calculate p-values using the null JPA-scores.
3) args.rr
    Performs RR-score analysis and calculate p-values from an analytical null model.
    Two analytical null models are implemented:
        a) --null maf
            MAF null (old deprecated version, not discussed in the manuscript)
            Requires 1000G / population MAF file, see documentation.
        b) --null perm
            PERM null (null model obtained by permuting SNPs)
'''

import numpy as np
import time
import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI

from utils.args import Args
from utils.logs import MyLogger
from qstats.jpa import JPA
from qstats.jpa_null import JPANULL
from qstats.revreg import RevReg

from iotools.data import Data
from iotools.outhandler import Outhandler
from iotools import readmaf
from iotools import readqnull


# ==================================
# Start MPI calculation
# ==================================
MPI.Init()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncore = comm.Get_size()
if rank == 0: start_time = time.time()


# ==================================
# Input Processing
# ==================================
args = Args(comm, rank)
logger = MyLogger(__name__)

# List of variables that are broadcast over all slave nodes
gtcent = None   # Centered genotype (x - mu), not divided by standard deviation, I x N
gtnorm = None   # Centered and scaled genotype (x - mu) / sig, dimension I x N
expr   = None   # Expression matrix, dimension G x N
masklist = None # List of gene indices masked for each SNP
maskcomp = None # List of CisMask (see utils/containers)
maf = None      # List of MAF of each SNP as observed in the sample (or read separately from the population if file is provided)
if rank == 0:
    logger.debug("Using {:d} cores".format(ncore))
    data = Data(args)
    data.load()
    gtcent = data.geno_centered
    gtnorm = data.geno_normed
    expr = data.expression
    masklist = data.cismasks_list
    maskcomp = data.cismasks_comp
    snpinfo = data.snpinfo # slaves don't need this. Only needed for maf + outhandler in master node
    geneinfo = data.geneinfo # slaves don't need this. Only needed for outhandler in master node
    maf = readmaf.load(snpinfo, maf_file = args.maf_file)
    logger.debug("After prefilter: {:d} SNPs and {:d} genes in {:d} samples".format(gtcent.shape[0], expr.shape[0], gtcent.shape[1]))
gtcent = comm.bcast(gtcent, root = 0)
gtnorm = comm.bcast(gtnorm, root = 0)
expr   = comm.bcast(expr,  root = 0)
masklist = comm.bcast(masklist, root = 0)
maskcomp = comm.bcast(maskcomp, root = 0)
maf  = comm.bcast(maf, root = 0)
comm.barrier()
if rank == 0: read_time = time.time()


# ==================================
# Calculate the null JPA-scores
# ==================================
if args.jpa_calc_null:
    if rank == 0: logger.debug("Calculating null values for JPA")
    jpa_null = JPANULL(gtnorm, expr, comm, rank, ncore, args.jpanull_file, niter=1000)
    jpa_null.compute()
if rank == 0: jpanull_time = time.time()


# ==================================
# Calculate JPA-score
# with empirical p-values from the null file
# TO-DO: 
#   - use covariate corrected gene expression
#   - why calculate full JPA-score for target genes of trans-eQTL SNPs?
# ==================================
if rank == 0: logger.debug("Computing JPA")
if args.jpa:
    jpa = JPA(gtnorm, expr, comm, rank, ncore, args.jpa, masklist, get_pvals = True, qnull = args.jpanull_file)
    jpa.compute()
if args.rr:
    # Run JPA for the target genes.
    jpa = JPA(gtnorm, expr, comm, rank, ncore, args.jpa, masklist, get_pvals = False)
    jpa.compute()
if rank == 0: jpa_time = time.time()


# ==================================
# Calculate the RR-score
# either with perm-null (default)
# or with MAF-null (old and deprecated, kept for legacy)
# TO-DO:
#  - The theory of MAF-null is not included in the manuscript, create a document
# ==================================
if args.rr:
    sigbeta2 = np.repeat(args.sigmabeta ** 2, gtnorm.shape[0])
    if rank == 0: logger.debug("Sigma_beta fixed to {:g}".format(args.sigmabeta))
    if args.nullmodel == 'maf':
        rr = RevReg(gtnorm, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel, maf = maf, masks = maskcomp)
    elif args.nullmodel == 'perm':
        rr = RevReg(gtcent, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel, maf = maf, masks = maskcomp)
    # Set get_betas = True to obtain the coefficients of multiple regression (stored in rr.betas)
    rr.compute(get_betas = False)
if rank == 0: rr_time = time.time()


# ==================================
# Write output
# ==================================
if rank == 0: 
    ohandle = Outhandler(args, snpinfo, geneinfo)
    if args.jpa:
        ohandle.write_jpa_out(jpa)
    if args.rr:
        ohandle.write_rr_out(jpa, rr)
if rank == 0: write_time = time.time()


# ==================================
# Log the run times
# ==================================
if rank == 0:
    logger.info("File reading time: {:g} seconds".format(read_time - start_time))
    logger.info("JPA null calculation time: {:g} seconds".format (jpanull_time - read_time))
    logger.info("JPA calculation time: {:g} seconds".format (jpa_time - jpanull_time))
    logger.info("RR calculation time: {:g} seconds".format (rr_time - jpa_time))
    logger.info("Result writing time: {:g} seconds".format(write_time - rr_time))
    logger.info("Total run time: {:g} seconds".format(time.time() - start_time))

MPI.Finalize()
